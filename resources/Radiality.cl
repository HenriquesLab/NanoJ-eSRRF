//#pragma OPENCL EXTENSION cl_khr_fp64: enable

float cubic(float x) {
    float a = 0.5f; // Catmull-Rom interpolation
    if (x < 0.0f) x = -x;
    float z = 0.0f;
    if (x < 1.0f)
        z = x * x * (x * (-a + 2.0f) + (a - 3.0f)) + 1.0f;
    else if (x < 2.0f)
        z = -a * x * x * x + 5.0f * a * x * x - 8.0f * a * x + 4.0f * a;
    return z;
}

float getInterpolatedValue(__global float* array, int const width, int const height, float const x, float const y) {
    int u0 = (int) floor(x - 0.5f);
    int v0 = (int) floor(y - 0.5f);
    float q = 0.0f;
    for (int j = 0; j <= 3; j++) {
        int v = min(max(v0 - 1 + j, 0), height-1);
        float p = 0.0f;
        for (int i = 0; i <= 3; i++) {
            int u = min(max(u0 - 1 + i, 0), width-1);
            p = p + array[v*width+u] * cubic(x - (u + 0.5f));
        }
        q = q + p * cubic(y - (v + 0.5f));
    }
    return q;
}

__kernel void calculateRadiality(
    __global float* xRingCoordinates,
    __global float* yRingCoordinates,
    __global float* pixels,
    __global float* GxArray,
    __global float* GyArray,
    __global float* radiality,
    __global float* interpolatedIntensity,
    int const magnification,
    float const spatialRadius,
    float const radialitySensitivity,
    float const shiftX,
    float const shiftY,
    int const nVectors,
    int const width,
    int const height,
    int const widthM,
    int const heightM
    )
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int offset = y * widthM + x;

    float xc = x + 0.5 + shiftX * magnification;
    float yc = y + 0.5 + shiftY * magnification;

    float CGH = 0; // CGH stands for Radiality original name - Culley-Gustafsson-Henriques transform

    float vx, vy, Gx, Gy;

    for (int vectorIter = 0; vectorIter < nVectors; vectorIter++) {
        vx = xc + xRingCoordinates[vectorIter];
        vy = yc + yRingCoordinates[vectorIter];

        Gx = getInterpolatedValue(GxArray, width, height, vx/magnification, vy/magnification);
        Gy = getInterpolatedValue(GyArray, width, height, vx/magnification, vy/magnification);
        float GMag = sqrt(Gx * Gx + Gy * Gy);

        // Calculate perpendicular distance from (xc,yc) to gradient line through (vx,vy)
        float absoluteDistance = fabs(Gy * (xc - vx) - Gx * (yc - vy)) / GMag;
        if (isnan(absoluteDistance)) absoluteDistance = spatialRadius;

        float Dk = 1 - absoluteDistance / spatialRadius; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx)
        Dk = pow(Dk, radialitySensitivity);

        // Accumulate Variables
        float GdotR = (Gx * xRingCoordinates[vectorIter] + Gy * yRingCoordinates[vectorIter]) / (GMag * spatialRadius); // tells you if vector was pointing inward or outward
        if (GdotR > 0 || GdotR != GdotR) CGH -= Dk; // vector was pointing outward
        else CGH += Dk; // vector was pointing inward
    }

//    if (CGH <= 0) { // less shadow, less funkiness
//        radiality[offset] = 0;
//        return;
//    }

    interpolatedIntensity[offset] = getInterpolatedValue(pixels, width, height, ((float) x)/magnification + shiftX, ((float) y)/magnification + shiftY);
    radiality[offset] = CGH / nVectors;
}


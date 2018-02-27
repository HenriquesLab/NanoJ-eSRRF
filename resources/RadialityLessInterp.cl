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

float getVBoundaryCheck(__global float* array, int const width, int const height, int const x, int const y) {
    int _x = min(max(x, 0), width-1);
    int _y = min(max(y, 0), height-1);
    return array[_y*width+_x];
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
    float distanceWeightSum = 0;

    float vx, vy, Gx, Gy;
    //float maxDistance = 1.4142135623730951;
    //float maxDistance = 2;

    for (int j=-1; j<=1; j++) {
        for (int i=-1; i<=1; i++) {
            vx = (int) (xc / magnification) + i + 0.5;
            vy = (int) (yc / magnification) + j + 0.5;

            Gx = getVBoundaryCheck(GxArray, width, height, vx, vy);
            Gy = getVBoundaryCheck(GyArray, width, height, vx, vy);
            float GMag = sqrt(Gx * Gx + Gy * Gy);
            //float distanceWeight = 1 - sqrt(pow(vx * magnification - xc, 2)+pow(vy * magnification - yc, 2)) / maxDistance;
            //float distanceWeight = 1/sqrt(pow(vx * magnification - xc, 2)+pow(vy * magnification - yc, 2));
            float distanceWeight = 1;

            // Calculate perpendicular distance from (xc,yc) to gradient line through (vx,vy)
            float absoluteDistance = fabs(Gy * (xc - vx * magnification) - Gx * (yc - vy * magnification)) / GMag;
            if (isnan(absoluteDistance)) absoluteDistance = spatialRadius;

            float Dk = 1 - absoluteDistance / spatialRadius; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx)
            Dk = pow(Dk, radialitySensitivity) * distanceWeight;
            distanceWeightSum += distanceWeight;

            // Accumulate Variables
            float GdotR = (Gx * i * magnification + Gy * j * magnification) / (GMag * spatialRadius); // tells you if vector was pointing inward or outward
            if (GdotR > 0 || GdotR != GdotR) CGH -= Dk; // vector was pointing outward
            else CGH += Dk; // vector was pointing inward
        }
    }

    interpolatedIntensity[offset] = getInterpolatedValue(pixels, width, height, ((float) x)/magnification + shiftX, ((float) y)/magnification + shiftY);
    radiality[offset] = CGH / distanceWeightSum;
}


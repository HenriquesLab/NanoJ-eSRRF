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
    int xM = get_global_id(0);
    int yM = get_global_id(1);
    int offset = yM * widthM + xM;

    float xc = (xM + 0.5) / magnification + shiftX;
    float yc = (yM + 0.5) / magnification + shiftY;

    float CGH = 0; // CGH stands for Radiality original name - Culley-Gustafsson-Henriques transform
    float distanceWeightSum = 0;

    float vx, vy, Gx, Gy;

    for (int j=-2; j<=2; j++) {
        for (int i=-2; i<=2; i++) {
            vx = (int) xc  + i + 0.5;
            vy = (int) yc  + j + 0.5;
            float maxDistance = fmax(sqrt(pow(vx - xc, 2)+pow(vy - yc, 2)),1);

            Gx = getVBoundaryCheck(GxArray, width, height, vx, vy);
            Gy = getVBoundaryCheck(GyArray, width, height, vx, vy);
            float GMag = sqrt(Gx * Gx + Gy * Gy);

            float distanceWeight = 1/maxDistance;

            // Calculate perpendicular distance from (xc,yc) to gradient line through (vx,vy)
            float absoluteDistance = fabs(Gy * (xc - vx) - Gx * (yc - vy)) / GMag;
            if (isnan(absoluteDistance)) absoluteDistance = maxDistance;

            float Dk = 1 - absoluteDistance / maxDistance; // Dk is now between 0 to 1, 1 if vector points precisely to (xc, yx)
            Dk *= distanceWeight;
            distanceWeightSum += distanceWeight;

            // Accumulate Variables
            float GdotR = (Gx * i * magnification + Gy * j * magnification); // tells you if vector was pointing inward or outward
            if (GdotR <= 0) CGH += Dk; // vector was pointing inwards
            else CGH -= Dk; // vector was pointing outwards
        }
    }
    CGH /= distanceWeightSum;
    if (CGH >= 0) CGH = fabs(pow(CGH, radialitySensitivity));
    else CGH = 0;

    interpolatedIntensity[offset] = getInterpolatedValue(pixels, width, height, ((float) xM)/magnification + shiftX, ((float) yM)/magnification + shiftY);
    radiality[offset] = CGH;
}


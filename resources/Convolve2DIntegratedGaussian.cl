#define ROOT2 1.4142135623730951f // square root of 2

// based on Smith et al, Nath Meth, 2010: Fast, single-molecule localization that achieves theoretically
// minimum uncertainty (see Sup Mat page 10)
// note, the paper has an error on their formula 4a and 4b, 2sigma^2 should be sqrt(2)*sigma
// see https://en.wikipedia.org/wiki/Normal_distribution formula 'Cumulative distribution function'
static float getIntegratedGaussian(float dx, float dy, float sigma2) {
    float Ex = 0.5f * (erf((dx + 0.5f) / sigma2) - erf((dx - 0.5f) / sigma2));
    float Ey = 0.5f * (erf((dy + 0.5f) / sigma2) - erf((dy - 0.5f) / sigma2));
    float vKernel = Ex * Ey;
    return vKernel;
}

__kernel void convolveHorizontal(
    __global float* pixels,
    __global float* pixelsConvolvedH,
    int const frameStart,
    int const frameStop,
    float const sigma) {

    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int w = get_global_size(0);
    const int h = get_global_size(1);
    const int wh = w * h;

    const int radius = max(((int) sigma) * 3, 1);
    const float sigma2 = ROOT2*fabs(sigma); // 2 * pow(sigma, 2);

    for (int f=frameStart; f<frameStop; f++) {
        int fOffset = f * wh;
        float vKernelSum = 0;
        float v = 0;

        for (int dx = -radius; dx <= radius; dx++) {
            float vKernel = 0.5f * (erf((dx + 0.5f) / sigma2) - erf((dx - 0.5f) / sigma2));
            int xp = min(max(x+dx, 0), w-1);
            v += pixels[fOffset + y * w + xp] * vKernel;
            vKernelSum += vKernel;
        }

        v /= vKernelSum;
        pixelsConvolvedH[fOffset + y * w + x] = v;
    }
}

__kernel void convolveVertical(
    __global float* pixelsConvolvedH,
    __global float* pixelsConvolved,
    int const frameStart,
    int const frameStop,
    float const sigma) {

    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int w = get_global_size(0);
    const int h = get_global_size(1);
    const int wh = w * h;

    const int radius = max(((int) sigma) * 3, 1);
    const float sigma2 = ROOT2*fabs(sigma); // 2 * pow(sigma, 2);

    for (int f=frameStart; f<frameStop; f++) {
        int fOffset = f * wh;
        float vKernelSum = 0;
        float v = 0;

        for (int dy = -radius; dy <= radius; dy++) {
            float vKernel = 0.5f * (erf((dy + 0.5f) / sigma2) - erf((dy - 0.5f) / sigma2));
            int yp = min(max(y+dy, 0), h-1);
            v += pixelsConvolvedH[fOffset + yp * w + x] * vKernel;
            vKernelSum += vKernel;
        }

        v /= vKernelSum;
        pixelsConvolved[fOffset + y * w + x] = v;
    }
}


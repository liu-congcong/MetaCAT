#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define white_x 95.047
#define white_y 100.0
#define white_z 108.883
#define gamma 2.4
#define eps (216.0 / 24389.0)
#define kappa (24389.0 / 27.0)
#define white_1x_15y_3z (white_x + 15.0 * white_y + 3.0 * white_z)
#define u (4.0 * white_x / white_1x_15y_3z)
#define v (9.0 * white_y / white_1x_15y_3z)
#define pi 3.14159265358979323846


static int calculate_xyz(double *xyz, double luminance, double chroma, double hue)
{
    /* xyz = x, y, z */
    *(xyz + 1) = luminance > eps * kappa ? pow((luminance + 16.0) / 116.0, 3.0) : luminance / kappa;
    double A = 1.0 / 3.0 * (52.0 * luminance / (chroma * cos(hue) + 13.0 * luminance * u) - 1.0);
    double B = -5.0 * (*(xyz + 1));
    *xyz = ((*(xyz + 1)) * (39.0 * luminance / (chroma * sin(hue) + 13.0 * luminance * v) - 5.0) - B) / (A + 1.0 / 3.0);
    *(xyz + 2) = (*xyz) * A + B;
    return 0;
}

static int xyz2rgb(double *xyz, double *rgb)
{
    double x;
    x = 3.2404542 * (*xyz) - 1.5371385 * (*(xyz + 1)) - 0.4985314 * (*(xyz + 2));
    rgb[0] = fmin(fmax(x > 0.0031308 ? 1.055 * pow(x, (1.0 / gamma)) - 0.055 : 12.92 * x, 0), 1);
    x = -0.9692660 * (*xyz) + 1.8760108 * (*(xyz + 1)) + 0.0415560 * (*(xyz + 2));
    rgb[1] = fmin(fmax(x > 0.0031308 ? 1.055 * pow(x, (1.0 / gamma)) - 0.055 : 12.92 * x, 0), 1);
    x = 0.0556434 * (*xyz) - 0.2040259 * (*(xyz + 1)) + 1.0572252 * (*(xyz + 2));
    rgb[2] = fmin(fmax(x > 0.0031308 ? 1.055 * pow(x, (1.0 / gamma)) - 0.055 : 12.92 * x, 0), 1);
    return 0;
}

double (*generateColors(int n))[3]
{
    double min_hue = 15.0;
    double max_hue = 375.0;
    double chroma = 100.0;
    double luminance = 65.0;
    double (*rgb)[3] = malloc(n * sizeof(double[3]));
    double xyz[3];
    if (!fmod(max_hue - min_hue, 360)) max_hue -= 360.0 / n;
    double hue_step = (max_hue - min_hue) / (n > 1 ? (n - 1) : 1);
    for (int i = 0; i < n; i++)
    {
        double hue = fmin(fmax(min_hue + i * hue_step, 0.0), 360.0) * pi / 180.0;
        calculate_xyz(xyz, luminance, chroma, hue);
        xyz2rgb(xyz, rgb[i]);
    }
    return rgb;
}

int freeColors(double (*x)[3])
{
    free(x);
    return 0;
}

// int printHelp()
// {
//     puts("hclColorGenerator v1.0.0");
//     puts("Generate colors in the HCL color space.");
//     puts("https://github.com/liu-congcong/hclColorGenerator");
//     puts("\nUsage:");
//     puts("  hclColorGenerator [options]");
//     puts("\nOptions:");
//     puts("  -n            Number of colors to generate");
//     puts("  -min-hue      Minimum hue (degrees, default: 15.0)");
//     puts("  -max-hue      Maximum hue (degrees, default: 375.0)");
//     puts("  -chroma       Chroma value (default: 100.0)");
//     puts("  -luminance    Luminance value (default: 65.0)");
//     return 0;
// }

// int main(int argc, char *argv[])
// {
//     int n = 0;
//     double min_hue = 15.0;
//     double max_hue = 375.0;
//     double chroma = 100.0;
//     double luminance = 65.0;

//     for (int i = 1; i < argc; i++)
//     {
//         if (!strncmp(argv[i], "-n", 2)) n = atoi(argv[i + 1]);
//         else if (!strncmp(argv[i], "-min-hue", 3)) min_hue = atof(argv[i + 1]);
//         else if (!strncmp(argv[i], "-max-hue", 3)) max_hue = atof(argv[i + 1]);
//         else if (!strncmp(argv[i], "-chroma", 2)) chroma = atof(argv[i + 1]);
//         else if (!strncmp(argv[i], "-luminance", 2)) luminance = atof(argv[i + 1]);
//     }
//     if (n <= 0)
//     {
//         printHelp();
//         exit(EXIT_FAILURE);
//     }

//     double (*rgb)[3] = malloc(n * sizeof(double[3]));
//     double xyz[3];
//     if (!fmod(max_hue - min_hue, 360)) max_hue -= 360.0 / n;
//     double hue_step = (max_hue - min_hue) / (n > 1 ? (n - 1) : 1);
//     for (int i = 0; i < n; i++)
//     {
//         double hue = fmin(fmax(min_hue + i * hue_step, 0.0), 360.0) * pi / 180.0;
//         calculate_xyz(xyz, luminance, chroma, hue);
//         xyz2rgb(xyz, rgb[i]);
//     }
//     for (int i = 0; i < n; i++) printf("%lf %lf %lf\n", rgb[i][0], rgb[i][1], rgb[i][2]);
//     freeColors(rgb);
//     return 0;
// }
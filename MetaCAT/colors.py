from collections.abc import Sequence
from math import cos, pi, sin

import numpy
from matplotlib import colormaps


def isValidColor(color):
    assert color in set(colormaps), f'\"{color}\" is not a valid color string.'
    return None


def parseRGBA(color: str, n: int, minColor: float, maxColor: float) -> Sequence[tuple[float, float, float, float]]:
    x = numpy.linspace(minColor, maxColor, n, endpoint = True)
    return colormaps[color](x)


def hcl(n, minHue = 15, maxHue = 375, chroma = 100, luminance = 65):
    colors = list()
    whiteX = 95.047
    whiteY = 100
    whiteZ = 108.883
    gamma = 2.4
    eps = 216 / 24389
    kappa = 24389 / 27
    xyz2rgb = ((3.2404542, -1.5371385, -0.4985314), (-0.9692660, 1.8760108, 0.0415560), (0.0556434, -0.2040259, 1.0572252))
    white1x15y3z = whiteX + 15 * whiteY + 3 * whiteZ
    u = 4 * whiteX / white1x15y3z
    v = 9 * whiteY / white1x15y3z
    y = pow((luminance + 16) / 116, 3) if (luminance > eps * kappa) else luminance / kappa
    b = - 5 * y
    if not (maxHue - minHue) % 360:
        maxHue -= 360 / n
    hueStep = (maxHue - minHue) / (n - 1) if n > 1 else 0
    for i in range(n):
        hue = min(max(minHue + i * hueStep, 0), 360) * pi / 180
        a = 1 / 3 * (52 * luminance / (chroma * cos(hue) + 13 * luminance * u) - 1)
        x = (y * (39 * luminance / (chroma * sin(hue) + 13 * luminance * v) - 5) - b) / (a + 1 / 3)
        z = x * a + b
        R = x * xyz2rgb[0][0] + y * xyz2rgb[0][1] + z * xyz2rgb[0][2]
        R = min(max(round(255.0 * (1.055 * pow(R, (1 / gamma)) - 0.055 if R > 0.0031308 else 12.92 * R)), 0), 255)
        G = x * xyz2rgb[1][0] + y * xyz2rgb[1][1] + z * xyz2rgb[1][2]
        G = min(max(round(255.0 * (1.055 * pow(G, (1 / gamma)) - 0.055 if G > 0.0031308 else 12.92 * G)), 0), 255)
        B = x * xyz2rgb[2][0] + y * xyz2rgb[2][1] + z * xyz2rgb[2][2]
        B = min(max(round(255.0 * (1.055 * pow(B, (1 / gamma)) - 0.055 if B > 0.0031308 else 12.92 * B)), 0), 255)
        colors.append(f'#{R:02X}{G:02X}{B:02X}')
    return colors

#ifndef FOURIER_H
#define FOURIER_H

// Calculate
// \int_{lower}^{upper} sin(w * x) dx
//   / cos(w * lower) / w - cos(w * upper) / w      if w != 0
// = |
    //   \ 0                                            if w == 0
inline double integral_sin(double w, double lower, double upper)
{
    return fabs(w) < 1e-16 ?
    0.0 :
        (cos(w * lower) - cos(w * upper)) / w;
}

// Calculate
// \int_{lower}^{upper} cos(w * x) dx
//   / sin(w * upper) / w - sin(w * lower) / w      if w != 0
// = |
    //   \ upper - lower                                if w == 0
inline double integral_cos(double w, double lower, double upper)
{
    return fabs(w) < 1e-16 ?
    upper - lower :
        (sin(w * upper) - sin(w * lower)) / w;
}

// Calculate
// \int_{lower}^{upper} sin(w1 * x) * cos(w2 * x) dx
// = 0.5 * \int sin((w1 - w2) * x) + sin((w1 + w2) * x) dx
inline double integral_sincos(double w1, double w2,
                              double lower, double upper)
{
    return 0.5 * (integral_sin(w1 - w2, lower, upper) +
                      integral_sin(w1 + w2, lower, upper));
}

// Calculate
// \int_{lower}^{upper} sin(w1 * x) * sin(w2 * x) dx
// = 0.5 * \int cos((w1 - w2) * x) - cos((w1 + w2) * x) dx
inline double integral_sinsin(double w1, double w2,
                              double lower, double upper)
{
    return 0.5 * (integral_cos(w1 - w2, lower, upper) -
                      integral_cos(w1 + w2, lower, upper));
}

// Calculate
// \int_{lower}^{upper} cos(w1 * x) * cos(w2 * x) dx
// = 0.5 * \int cos((w1 - w2) * x) + cos((w1 + w2) * x) dx
inline double integral_coscos(double w1, double w2,
                              double lower, double upper)
{
    return 0.5 * (integral_cos(w1 - w2, lower, upper) +
                      integral_cos(w1 + w2, lower, upper));
}

#endif  // FOURIER_H

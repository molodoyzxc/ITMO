#include <math.h>
#include <stdio.h>

struct Electron {
    double x, y, r, R, dx, dy, q, m, L, t;
};

void InitElectron(struct Electron *electron, double r, double R, double dx, double L) {
    electron->x = 0;
    electron->r = r;
    electron->R = R;
    electron->dx = dx;
    electron->dy = 0;
    electron->L = L;
    electron->y = (R - r) / 2 + r;
    electron->q = -1.6 * pow(10, -19);
    electron->m = 9.1 * pow(10, -31);
    electron->t = 0;
}

void Motion(struct Electron *electron, double U) {
    double dt = 0.000000000001;
    while (electron->x < electron->L & electron->y > electron->r) {
        double a = (electron->q * U) / (electron->y * electron->m * log(electron->R / electron->r));
        electron->dy += a * dt;
        electron->y += electron->dy * dt;
        electron->x += electron->dx * dt;
        electron->t += dt;
    }
}

void MotionGraphic(struct Electron *electron, double U) {
    double dt = 0.000000000001;
    FILE *y, *a, *t, *x, *v;
    y = fopen("graphs\\y.txt", "w");
    v = fopen("graphs\\v.txt", "w");
    a = fopen("graphs\\a.txt", "w");
    t = fopen("graphs\\t.txt", "w");
    x = fopen("graphs\\x.txt", "w");
    while (electron->x < electron->L | electron->y > electron->r) {
        fprintf(y, "%f\n", electron->y);
        fprintf(x, "%f\n", electron->x);
        fprintf(t, "%f\n", electron->t * 1000000000);
        fprintf(v, "%f\n", electron->dy);
        double ddy = (electron->q * U) / (electron->y * electron->m * log(electron->R / electron->r));
        fprintf(a, "%f\n", ddy);
        electron->dy += ddy * dt;
        electron->y += electron->dy * dt;
        electron->x += electron->dx * dt;
        electron->t += dt;
    }
    fclose(x);
    fclose(y);
    fclose(a);
    fclose(t);
    fclose(v);
}

void ChangeComa(FILE *in, FILE *out) {
    int c;
    while ((c = fgetc(in)) != EOF) {
        if (c == '.') {
            c = ',';
        }
        fputc(c, out);
    }
    fclose(in);
    fclose(out);
}

void Convert() {
    FILE *in, *out;
    in = fopen("graphs\\x.txt", "r");
    out = fopen("graphs\\x.txt", "r+");
    ChangeComa(in, out);
    in = fopen("graphs\\y.txt", "r");
    out = fopen("graphs\\y.txt", "r+");
    ChangeComa(in, out);
    in = fopen("graphs\\t.txt", "r");
    out = fopen("graphs\\t.txt", "r+");
    ChangeComa(in, out);
    in = fopen("graphs\\a.txt", "r");
    out = fopen("graphs\\a.txt", "r+");
    ChangeComa(in, out);
    in = fopen("graphs\\v.txt", "r");
    out = fopen("graphs\\v.txt", "r+");
    ChangeComa(in, out);
}

int main() {
    double Umax = 1000, Umin = 0, U;
    struct Electron electron;
    while (Umax - Umin > 0.0000001) {
        InitElectron(&electron, 0.1, 0.21, 850000, 0.29);
        U = (Umax + Umin) / 2;
        Motion(&electron, U);
        if (electron.x >= electron.L) {
            Umin = U;
        } else {
            Umax = U;
        }
    }
    InitElectron(&electron, 0.1, 0.21, 850000, 0.29);
    MotionGraphic(&electron, U);
    printf("U %f V\nt %f nanosec\nV %f m/s\n", U, electron.t * 1000000000,
           sqrt(pow(electron.dy, 2) + pow(electron.dx, 2)));
    return 0;
}
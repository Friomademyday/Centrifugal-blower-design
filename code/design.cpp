#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <cstdio>
#include <conio.h>

const float PI = 3.1429f;
const float RAD_TO_DEG = 180.0f / PI;
const float DEG_TO_RAD = PI / 180.0f;
const float G = 9.81f;

float Delta_p;
float Q;
float n;
float Had;
float sigma, sigma1;
float Ov_Eff, Hy_Eff, Vo_Eff, Me_Eff;
float N_o, N_i, N_imax;
float Qth;
float C0;
float dh;
float D0;
float Km0;
float Cm0;
float A;
float D1;
float u1;
float Km1;
float Cm1;
float Beta1;
float dlt;
float Z;
float t1;
float Km11;
float Km12;
float DeltaBeta;
float Beta1f;
float t1f;
float Km1f;
float B1;
float Mu;
float B11;
float X;
float Cm2;
float w;
float temp1;
float Beta2;
float E;
float C_factor_cu2_u2;
float u2;
float D2;
float D1_dash;
float temp2;
float Y;
float Si1;
float p;
float Delta_pTinf;
float F_factor;
float u21;
float D21;
float D1_dash1;
float temp3;
float p1;
float Delta_pTinf1;
float F1_factor;
float u22;
float D22;
float D1_dash2;
float temp4;
float p2;
float Delta_pTinf2;
float F2_factor;
float u23;
float D23;
float t2;
float Km2;
float B2;
float Cm3;
float Zcal;
float gamma_val = 1.225f;
float rho = 0.125f;

float Sdr_imp = 0.0f;
float dr_imp = 0.0f;
float w1_imp;
float w2_imp;
float Cm_imp;
float tt_imp;
float t_imp;
float Beta_imp;
float F_imp;
float Cmr_imp;
float B_imp;
float Sd0_imp = 0.0f;
float d0_imp = 0.0f;
float AvgG_imp;

float Cu23_vol;
float Cu3_vol;
float alpha3_vol;
float C3_vol;
float p_vol;
float B3_vol;
float R3_vol;

std::vector<float> r_imp;
std::vector<float> G_imp;
std::vector<int> theta_vol = {0, 45, 90, 135, 180, 225, 270, 315, 360};
std::vector<float> Q_theta_vol;
std::vector<float> area_theta_vol;
std::vector<float> height_theta_vol;
std::vector<float> height_eff_vol;

void print_heading();
void head_impeller();
void head_impl();
void impeller();
void head_volute();
void volute();

int main() {
    printf("\n\n\t\t DESIGN OF CENTRIFUGAL BLOWER - \n");

    Delta_p = 140.0f;
    Q = 1.61f;
    n = 2880.0f;
    Had = Delta_p / gamma_val;

    sigma = ((0.379f * n * sqrt(Q)) / (60.0f * pow(Had, 0.75f)));
    sigma1 = ((n * sqrt(Q)) / (28.5f * pow((Delta_p / rho), 0.75f)));

    printf("\n\n\t SPECIFICATIONS ");
    printf("\n\n Delta_p = %3.0f mmwc", Delta_p);
    printf("\n\n Q = %2.2f m3/sec ", Q);
    printf("\n\nn = %d rpm", static_cast<int>(n));
    printf("\n\n Adiabatic Head (Had) = %2.1f m", Had);
    printf("\n\n speed coefficient(sigma) = %1.2f", sigma);
    printf("\n\n speed coefficient(sigma1) = %1.2f", sigma1);

    printf("\n\n\tEnter the value of Me_Eff (0.92 to 0.98):");
    scanf("%f", &Me_Eff);
    printf("\n\n\tEnter the value of Hy_Eff (0.76 to 0.86):");
    scanf("%f", &Hy_Eff);
    printf("\n\n\tEnter the value of Vo_Eff (0.95 to 0.98):");
    scanf("%f", &Vo_Eff);

    Ov_Eff = Hy_Eff * Vo_Eff * Me_Eff;
    printf("\n\n Overall Efficiency = %0.3f", Ov_Eff);
    getch();

    N_o = 9.81f * Delta_p * Q / 1000.0f;
    N_i = N_o / Ov_Eff;
    N_imax = 1.1f * N_i;

    printf("\n\n\t N_o: %0.2f kW", N_o);
    printf("\n\n\t N_i: %0.2f kW", N_i);
    printf("\n\n\t N_imax : %0.2f kW", N_imax);

    Qth = Q / Vo_Eff;

    printf("\n\n\tSUCTION CONDITIONS ");
    C0 = 0.095f * pow((Qth * n * n), 0.333f);

    printf("\n\n\tEnter the value of hub diameter dh in m:");
    scanf("%f", &dh);

    temp1 = (4.0f * Qth) / (PI * C0);
    D0 = sqrt(temp1 + (dh * dh));

    printf("\n\n\t CO %3.2f m/sec", C0);
    printf("\n\n\t D0 : %1.3f m", D0);

    printf("\n\n\t\t INLET PARAMETERS");
    printf("\n\n\tEnter the value of Km0(=Cm0/C0) (0.9 to 1.0):");
    scanf("%f", &Km0);

    Cm0 = Km0 * C0;
    printf("\n\n\t Cm0: %3.2f m/sec", Cm0);

    printf("\n\n\tEnter the value of A in m (0.000 to 0.002):");
    scanf("%f", &A);

    D1 = D0 + A;
    u1 = (PI * D1 * n) / 60.0f;

    printf("\n\n\t\t Ist APPROXIMATION");
    printf("\n\n\tEnter the value of Km1(=Cm1/Cm0) (1.0 to 1.5):");
    scanf("%f", &Km1);

    Cm1 = Km1 * Cm0;
    Beta1 = atan(Cm1 / u1);
    printf("\n\n\t Beta1 = %2.1f degrees", Beta1 * RAD_TO_DEG);

    printf("\n\n\t\t IInd APPROXIMATION");
    printf("\n\n\tEnter the value of dlt in m:");
    scanf("%f", &dlt);
    printf("\n\n\tEnter the value of Z :");
    scanf("%f", &Z);

    t1 = PI * D1 * sin(Beta1);
    Km11 = 1.0f / (1.0f - ((Z * dlt) / t1));
    Cm1 = Km11 * Cm0;
    Beta1 = atan(Cm1 / u1);
    printf("\n\n\t Beta1 = %2.1f degrees", Beta1 * RAD_TO_DEG);

    printf("\n\n\t\t IIIrd APPROXIMATION");
    t1 = PI * D1 * sin(Beta1);
    Km12 = 1.0f / (1.0f - ((Z * dlt) / t1));
    Cm1 = Km12 * Cm0;
    Beta1 = atan(Cm1 / u1);
    printf("\n\n\t Beta1 = %2.1f degrees", Beta1 * RAD_TO_DEG);

    printf("\n\n\tEnter the value of DeltaBeta in degrees:");
    scanf("%f", &DeltaBeta);

    DeltaBeta = (DeltaBeta * DEG_TO_RAD);
    Beta1f = Beta1 + DeltaBeta;
    t1f = PI * D1 * sin(Beta1f);
    Km1f = 1.0f / (1.0f - ((Z * dlt) / t1f));
    Cm1f = Km1f * Cm0;
    B1 = (Qth / (PI * D1 * Cm0));

    printf("\n\n\t Enter the value of Mu(0.76 to 0.97) :");
    scanf("%f", &Mu);

    B11 = D1 / (4.0f * Mu);

    printf("\n\n\t\t INLET VALUES ");
    printf("\n\n Qth =%1.2f m3/sec", Qth);
    printf("\n\n D0 =%1.3f m", D0);
    printf("\n\n D1 =%1.3f m", D1);
    printf("\n\n u1 =%2.2f m/sec", u1);
    printf("\n\n B1 =%1.3f m", B1);
    printf("\n\n B11 =%1.3f m", B11);
    printf("\n\n CO =%2.2f m/sec", C0);
    printf("\n\n Cm0 = %2.2f m/sec", Cm0);
    printf("\n\n Cm1 = %2.2f m/sec", Cm1);
    printf("\n\n Cm1f = %2.2f m/sec", Cm1f);
    printf("\n\n Beta1f=%2.1f degrees", Beta1f * RAD_TO_DEG);
    printf("\n\n Km1f = %1.2f", Km1f);
    getch();

    printf("\n\n\t OUTLET PARAMETERS\n");
    printf("\n\n\t Enter the value of(Cm_dash=Cm2/Cm1) X(0.8 to 0.9):");
    scanf("%f", &X);
    Cm2 = (X * Cm1f);
    printf("\n\n Cm2 =%2.2f m/sec", Cm2);

    printf("\n\n\t Enter the value of(w_dash=w1/w2) w(1.25 to 1.4):");
    scanf("%f", &w);
    temp1 = w * X * sin(Beta1f);
    Beta2 = asin(temp1);
    printf("\n\n Beta2 =%2.1f degrees", Beta2 * RAD_TO_DEG);

    E = Cm2 / (2.0f * tan(Beta2));
    printf("\n\n\t E = %2.2f ", E);

    printf("\n\n\t Enter the value of (Cu_dash=cu2/u2) C (0.5 to 0.8):");
    scanf("%f", &C_factor_cu2_u2);

    u2 = sqrt(Delta_p / (C_factor_cu2_u2 * rho * Hy_Eff));
    D2 = (60.0f * u2) / (PI * n);
    printf("\n\n u2 =%2.2f m/sec", u2);
    printf("\n\n D2 =%1.3f m", D2);

    printf("\n\n\t Ist APPROXIMATION\n");
    D1_dash = (D1 / D2);
    printf("\n\n D1_dash =%2.2f", D1_dash);
    temp2 = (D1_dash * D1_dash);

    printf("\n\n\t Enter the value of Y (0.5 to 0.65):");
    scanf("%f", &Y);

    Si1 = Y + 0.6f * sin(Beta2);
    p = (2.0f * Si1 / Z) * (1.0f / (1.0f - temp2));
    printf("\n\n Si1 =%1.2f", Si1);
    printf("\n\np =%1.2f", p);

    Delta_pTinf = ((1.0f + p) * Delta_p) / (Hy_Eff);
    printf("\n\n\t Delta_pTinf = %4.2f", Delta_pTinf);
    F_factor = Delta_pTinf / rho;
    printf("\n\n\t F = %4.2f", F_factor);

    u21 = E + sqrt((E * E) + F_factor);
    D21 = (60.0f * u21) / (PI * n);
    printf("\n\n u21 =%2.2f m/sec", u21);
    printf("\n\n D21 =%1.4f m", D21);

    printf("\n\n\t IInd APPROXIMATION\n");
    D1_dash1 = (D1 / D21);
    printf("\n\n D1_dash1 =%2.2f", D1_dash1);
    temp3 = (D1_dash1 * D1_dash1);

    p1 = (2.0f * Si1 / Z) * (1.0f / (1.0f - temp3));
    printf("\n\n p1 =%1.2f", p1);

    Delta_pTinf1 = ((1.0f + p1) * Delta_p) / (Hy_Eff);
    printf("\n\n\t Delta_pTinf1 = %4.2f ", Delta_pTinf1);
    F1_factor = Delta_pTinf1 / rho;
    printf("\n\n\t F1 = %4.2f", F1_factor);

    u22 = E + sqrt((E * E) + F1_factor);
    D22 = (60.0f * u22) / (PI * n);
    printf("\n\n u22 =%2.2f m/sec", u22);
    printf("\n\n D22 =%1.4f m", D22);

    printf("\n\n\t IIIrd APPROXIMATION\n");
    D1_dash2 = (D1 / D22);
    printf("\n\n D1_dash2 =%2.2f", D1_dash2);
    temp4 = (D1_dash2 * D1_dash2);

    p2 = (2.0f * Si1 / Z) * (1.0f / (1.0f - temp4));
    printf("\n\n p2 =%1.2f", p2);

    Delta_pTinf2 = ((1.0f + p2) * Delta_p) / (Hy_Eff);
    printf("\n\n\t Delta_pTinf2 = %4.2f", Delta_pTinf2);
    F2_factor = Delta_pTinf2 / rho;
    printf("\n\n\t F2 = %4.2f", F2_factor);

    u23 = E + sqrt((E * E) + F2_factor);
    D23 = (60.0f * u23) / (PI * n);
    printf("\n\n u23 =%2.2f m/sec", u23);
    printf("\n\n D23 =%1.4f m", D23);

    t2 = PI * D23 * sin(Beta2);
    Km2 = 1.0f / (1.0f - ((Z * dlt) / t2));
    B2 = Q / (PI * D23 * Cm2);
    Cm3 = Km2 * Cm2;
    Zcal = PI * 1.4f * ((1.0f + D1_dash2) / (1.0f - D1_dash2));

    printf("\n\n B2 =%1.3f m", B2);
    printf("\n\n\t\t Final Values ");
    printf("\n\nZ = %d ", static_cast<int>(Z));
    printf("\n\n Z calculated = %d ", static_cast<int>(Zcal));
    printf("\n\n B1 = %1.3f m", B1);
    printf("\n\n B2 = %1.3f m", B2);
    printf("\n\n C0 = %2.2f m/sec", C0);
    printf("\n\n Cm0 = %2.2f m/sec", Cm0);
    printf("\n\n Cm1f = %2.2f m/sec", Cm1f);
    printf("\n\n Cm2 = %2.2f m/sec", Cm2);
    printf("\n\n Cm3 =%1.2f m/sec ", Cm3);
    printf("\n\nD1 = %1.3f m ", D1);
    printf("\n\n D2 = %1.3f m ", D2);
    printf("\n\n D21 = %1.3f m", D21);
    printf("\n\n D22 = %1.3f m", D22);
    printf("\n\n D23 = %1.3f m", D23);
    printf("\n\n u2 = %2.1f m/sec ", u2);
    printf("\n\n u21 = %2.1f m/sec ", u21);
    printf("\n\n u22 = %2.1f m/sec ", u22);
    printf("\n\n u23 = %2.1f m/sec ", u23);
    printf("\n\n Beta1f = %2.1f degrees", Beta1f * RAD_TO_DEG);
    printf("\n\n Beta2 = %2.1f degrees", Beta2 * RAD_TO_DEG);

    impeller();
    volute();
    getch();
    return 0;
}

void impeller() {
    printf("\n\n VANE DEVELOPMENT -(SINGLE CURVATURE)\n");

    int nn_imp = 0;
    printf("\n\n\tEnter the no of values (7 to 10):");
    scanf("%d", &nn_imp);

    r_imp.resize(nn_imp);
    G_imp.resize(nn_imp);

    printf("\n\n\tEnter Radius r in m =\n\t");
    for (int i = 0; i < nn_imp; i++) {
        scanf("%f", &r_imp[i]);
    }
    getch();

    head_impeller();
    Sdr_imp = 0.0f;
    int count_imp = 1;

    if (nn_imp < 2) {
        printf("Error: Not enough radius values for impeller calculation.\n");
        return;
    }

    for (int i = 0; i < nn_imp; i++) {
        if (i < nn_imp - 1) {
            dr_imp = r_imp[i + 1] - r_imp[i];
            Sdr_imp += dr_imp;
        } else {
            dr_imp = 0.0f;
        }

        w1_imp = Cm1f / sin(Beta1f);
        w2_imp = Cm2 / sin(Beta2);
        Cm_imp = Cm1f + ((Cm2 - Cm1f) / ((D22 - D1) / 2.0f) * Sdr_imp);
        tt_imp = (w2_imp - w1_imp) / ((D22 - D1) / 2.0f);
        w_imp = w1_imp + (tt_imp * Sdr_imp);
        t_imp = (2.0f * PI * r_imp[i]) / Z;
        temp1 = Cm_imp / w_imp;

        if (std::abs(temp1) > 1.0f) {
            printf("\nWarning: arcsin argument out of range (Cm/w = %f). Clamping to +/-1.", temp1);
            temp1 = std::copysign(1.0f, temp1);
        }
        Beta_imp = asin(temp1);

        if (i == 0) {
            Beta_imp = Beta1f;
        }
        if ((i == 5) || (i == 6)) {
            Beta_imp = Beta2;
        }

        F_imp = tan(Beta_imp);
        if (r_imp[i] * F_imp == 0) {
            G_imp[i] = 0;
        } else {
            G_imp[i] = 1.0f / (r_imp[i] * F_imp);
        }

        Cmr_imp = Cm0 + ((Cm0 - Cm3) / ((D1 - D22) / 2.0f) * Sdr_imp);
        if (PI * r_imp[i] * Cmr_imp == 0) {
            B_imp = 0;
        } else {
            B_imp = Qth / (2.0f * PI * r_imp[i] * Cmr_imp);
        }

        printf("\n %d %1.4f %1.4f %1.2f %1.4f %2.2f %1.4f", count_imp, r_imp[i], Sdr_imp, Cm_imp, w_imp, (Beta_imp * RAD_TO_DEG), B_imp);
        count_imp++;
    }
    printf("\n-----------------------------------\n");
    getch();

    head_impl();
    Sd0_imp = 0.0f;
    d0_imp = 0.0f;

    for (int i = 0; i < nn_imp; i++) {
        if (i < nn_imp - 1) {
            AvgG_imp = (G_imp[i] + G_imp[i + 1]) / 2.0f;
            dr_imp = r_imp[i + 1] - r_imp[i];
            d0_imp = AvgG_imp * dr_imp;
        } else {
            AvgG_imp = G_imp[i];
            dr_imp = 0.0f;
            d0_imp = 0.0f;
        }

        printf("\n %1.4f %3.2f", G_imp[i], (Sd0_imp * RAD_TO_DEG));
        if (i < nn_imp - 1) {
            printf("\n\t\t\t %1.4f %1.4f %1.4f", AvgG_imp, dr_imp, d0_imp);
        }
        Sd0_imp += d0_imp;
    }
    printf("\n-----------------------------------\n");
    getch();
}

void volute() {
    printf("\n\n VOLUTE DESIGN\n");
    printf("\n\n CIRCULAR CROSS SECTION AND CONSTANT VELOCITY \n");

    Cu23_vol = ((9081.0f * F2_factor) / u23);
    Cu3_vol = Cu23_vol;
    temp1 = (Cm3 / Cu3_vol);

    if (std::abs(temp1) > 1.0f) {
        printf("\nWarning: arcsin argument out of range (Cm3/Cu3 = %f). Clamping to +/-1.", temp1);
        temp1 = std::copysign(1.0f, temp1);
    }
    alpha3_vol = atan(temp1);

    C3_vol = Cm3 / sin(alpha3_vol);

    printf("\n\n\tEnter the value of p (1.03 to 1.05):");
    scanf("%f", &p_vol);

    B3_vol = B2 + 0.05f * D23;
    R3_vol = p_vol * D23 / 2.0f;

    printf("\n\n R3 = %1.3f m", R3_vol);
    printf("\n\n B3 = %1.3f m", B3_vol);
    printf("\n\n alpha3 = %2.2f Degrees", alpha3_vol * RAD_TO_DEG);
    printf("\n\n C3 = %1.3f m/sec", C3_vol);
    getch();

    head_volute();

    Q_theta_vol.resize(theta_vol.size());
    area_theta_vol.resize(theta_vol.size());
    height_theta_vol.resize(theta_vol.size());
    height_eff_vol.resize(theta_vol.size());

    for (size_t i = 0; i < theta_vol.size(); i++) {
        Q_theta_vol[i] = ((Q / 360.0f) * theta_vol[i]);
        area_theta_vol[i] = Q_theta_vol[i] / Cm3;
        height_theta_vol[i] = area_theta_vol[i] / B3_vol;
        height_eff_vol[i] = 1.1f * height_theta_vol[i];

        printf("\n %zu %d %1.4f %1.4f %1.4f %1.4f \n", i + 1, theta_vol[i], Q_theta_vol[i], area_theta_vol[i], height_theta_vol[i], height_eff_vol[i]);
    }
    printf("\n-----------------------------------\n");
    getch();
}

void print_heading() {
    printf("-----------------------------------\n");
    printf("\n\n");
    printf(" Sr.No. Cm1 betal w1 K1 \n");
    printf("-----------------------------------\n");
}

void head_impeller() {
    printf("\n\n");
    printf("-----------------------------------\n");
    printf(" SNO     r        Sdr       Cm       Cmr      w        beta       B\n");
    printf("-----------------------------------\n");
}

void head_impl() {
    printf("\n\n");
    printf("-----------------------------------\n");
    printf(" G       Avg(G)   dr       d0        SdO \n");
    printf("-----------------------------------\n");
}

void head_volute() {
    printf("-----------------------------------\n");
    printf(" s.no_theta  theta  Q(theta) A(theta) h(theta) h(eff) \n");
    printf("-----------------------------------\n");
}

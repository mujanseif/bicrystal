const int n_nodes  = 4;
const int n_limits = 8;
const double pi = 4.0 * atan(1.0);

double init_point(double *i_vec1, double *i_vec2,
                  double *i_vec3, double *i_vec4,
                  int     i_vec_size){
  // Initialises points on the surface element given four vectors.
  double result = 0.0;
  double denom = 0.0;
  denom = dot_product(i_vec3, i_vec4, i_vec_size);
  if (denom == 0.0){
    fprintf(stderr, "nodal_surface_force_linear_rectangle: init_point: division by zero, dot_product(i_vec3, i_vec4) = %2.14f\n", denom);
    //exit(EXIT_FAILURE);
  }
  result = dot_product(i_vec1, i_vec2, i_vec_size)/denom;
  return result;
}


double seed_single_integral(double a, double b){
  // Single integral seed.
  double result = 0.0;
  double arg;
  arg = a+b;
  if (arg <= 0.0){
    fprintf(stderr, "nodal_surface_force_linear_rectangle: seed_single_integral: log(%2.14f) = %2.14f\n", arg, log(arg));
    //exit(EXIT_FAILURE);
  }
  result = log(arg);
  return result;
}

double integral_type_1(double a, double b,
                       double c, double d,
                       double e){
  double result = 0.0;
  result = 0.5*(a*b + (c-d)*e);
  return result;
}

double numer_seed_double_integral(double a,
                                  double b, double c, double d,
                                  double e){
  // Returns the argument of the numerator of the seed integral for double integrals.
  double result = 0.0;
  result = a*(b - c + d) + e;
  return result;
}

double denom_seed_double_integral(double a, double b,
                                  double c, double d){
  // Returns the argument of the denominator of the seed integral for double integrals.
  double result = 0.0;
  result = a*b + c*d;
  if(result == 0.0){
    fprintf(stderr, "nodal_surface_force_linear_rectangle: denom_seed_integral_00m3: division by 0, result = %2.14f\n", result);
    //exit(EXIT_FAILURE);
  }
return result;
}

double seed_double_integral(double a, double b, double c, double d, double e,
                            double f, double g, double h, double i){
  // Returns the seed integral for double integrals.
  double result      = 0.0;
  double numerator   = 0.0;
  double denominator = 0.0;
  numerator   = numer_seed_double_integral(a, b, c, d, e);
  denominator = denom_seed_double_integral(f, g, h, i);
  if(denominator > 0.){
    denominator = sqrt(denominator);
    result = 2.0/denominator * atan(numerator/denominator);
  }
  else{
    denominator = sqrt(-denominator);
    result = -2.0/denominator * atanh(numerator/denominator);
  }
  return result;
}

double integral_type_2(double a,
                       double b, double c){
  double result = 0.0;
  result = a + b*c;
  return result;
}

double integral_type_3(double a,
                       double b, double c,
                       double d, double e){
  double result = 0.0;
  result = a + b*c + d*e;
  return result;
}

double integral_type_4(double a,
                       double b, double c,
                       double d, double e,
                       double f, double g){
  double result = 0.0;
  result = a + b*c + d*e + f*g;
  return result;
}

double integral_type_5(double a,
                       double b, double c,
                       double d, double e,
                       double f, double g,
                       double h, double i){
  double result = 0.0;
  result = a + b*c + d*e + f*g + h*i;
  return result;
}

double integral_type_6(double a, double b,
                       double c, double d,
                       double e, double f, double g,
                       double h, double i,
                       double j, double k){
  double result = 0.0;
  result = a*b + c*d + (e+f)*g + h*i + j*k;
  return result;
}

double integral_type_7(double a, double b,
                       double c, double d,
                       double e, double f, double g,
                       double h, double i){
  double result = 0.0;
  result = a*b + c*d + (e+f)*g + h*i;
  return result;
}

double integral_type_8(double a, double b,
                       double c, double d,
                       double e, double f,
                       double g, double h){
  double result = 0.0;
  result = a*b + c*d + e*f + g*h;
  return result;
}

double integral_type_9(double a,
                       double b, double c,
                       double d, double e,
                       double f, double g,
                       double h, double i,
                       double j, double k){
  double result = 0.0;
  result = a + b*c + d*e + f*g + h*i + j*k;
  return result;
}

double integral_type_10(double a,
                        double b, double c,
                        double d, double e,
                        double f, double g,
                        double h, double i,
                        double j, double k){
  double result = 0.0;
  result = a + b*c + d*e + f*g + h*i + j*k;
  return result;
}

double integral_type_11(double a,
                        double b, double c,
                        double d, double e,
                        double f, double g,
                        double h, double i){
  double result = 0.0;
  result = a + b*c + d*e + f*g + h*i;
  return result;
}

double integral_type_12(double a, double b,
                        double c, double d,
                        double e, double f,
                        double g, double h,
                        double i, double j,
                        double k, double l, double m){
  double result = 0.0;
  result = a*b + c*d + e*f + g*h + i*j + k*l*m;
  return result;
}

void integral_vector(double *i_p, double *i_q, double *i_b, double *i_t, double *i_n,
                     double i_one_m_nu, double i_a_sq,
                     double o_vec_int[][3]){
  // Calculates the vectors that are multiplied by the integrals.
  // Has to be void in order to return a 2D array, otherwise we'd have to use dynamic memory allocation and therefore pointers. We don't need such flexibility as this is quite specific.
  /*
    (t x b), (p x b)
    (q x b), (b x t)
    t*((p x b) dot n), t*((q x b) dot n)
    t*((t x b) dot n), t*((b x t) dot n)
    (t dot n) * (p x b), (t dot n) * (q x b)
    (t dot n) * (t x b)
    n * ((p x b) dot t), n * ((q x b) dot t)
  */
  double t_x_b[3], p_x_b[3],
         q_x_b[3], b_x_t[3],
         t_p_x_b_dot_n[3], t_q_x_b_dot_n[3],
         t_t_x_b_dot_n[3], t_b_x_t_dot_n[3],
         t_dot_n_p_x_b[3], t_dot_n_q_x_b[3],
         t_dot_n_t_x_b[3],
         n_p_x_b_dot_t[3], n_q_x_b_dot_t[3];
  /*
    t dot n
    ((p x b) dot n), ((q x b) dot n)
    ((t x b) dot n), ((b x t) dot n)
    ((p x b) dot t), ((q x b) dot t)
    (t dot n) * ((p x b) dot t), (t dot n) * ((q x b) dot t)
    3*(t dot n) * ((p x b) dot t), 3*(t dot n) * ((q x b) dot t)
    1.5 * (1-nu) * a^2, a^3
  */
  double t_dot_n,
         p_x_b_dot_n, q_x_b_dot_n,
         t_x_b_dot_n, b_x_t_dot_n,
         p_x_b_dot_t, q_x_b_dot_t,
         t_dot_n_p_x_b_dot_t, t_dot_n_q_x_b_dot_t,
         t_dot_n_p_x_b_dot_t_3, t_dot_n_q_x_b_dot_t_3,
         one_m_nu_1p5_a_sq, a_sq_3;
  int i;
  // Calculate auxiliary local variables.
  one_m_nu_1p5_a_sq = 1.5 * i_one_m_nu * i_a_sq;
  a_sq_3 = 3.0 * i_a_sq;
  cross_product(i_p, i_b, p_x_b);
  cross_product(i_q, i_b, q_x_b);
  cross_product(i_t, i_b, t_x_b);
  cross_product(i_b, i_t, b_x_t);
  t_dot_n     = dot_product(i_t  , i_n, 3);
  p_x_b_dot_n = dot_product(p_x_b, i_n, 3);
  q_x_b_dot_n = dot_product(q_x_b, i_n, 3);
  t_x_b_dot_n = dot_product(t_x_b, i_n, 3);
  b_x_t_dot_n = dot_product(b_x_t, i_n, 3);
  p_x_b_dot_t = dot_product(p_x_b, i_t, 3);
  q_x_b_dot_t = dot_product(q_x_b, i_t, 3);
  t_dot_n_p_x_b_dot_t   = t_dot_n * p_x_b_dot_t;
  t_dot_n_q_x_b_dot_t   = t_dot_n * q_x_b_dot_t;
  t_dot_n_p_x_b_dot_t_3 = 3.0 * t_dot_n_p_x_b_dot_t;
  t_dot_n_q_x_b_dot_t_3 = 3.0 * t_dot_n_q_x_b_dot_t;

  for (i=0; i<3; i++)
  {
      // Preparing common array factors.
      t_t_x_b_dot_n[i] = i_t  [i]*t_x_b_dot_n;
      t_b_x_t_dot_n[i] = i_t  [i]*b_x_t_dot_n;
      t_p_x_b_dot_n[i] = i_t  [i]*p_x_b_dot_n;
      t_q_x_b_dot_n[i] = i_t  [i]*q_x_b_dot_n;
      t_dot_n_t_x_b[i] = t_x_b[i]*t_dot_n;
      t_dot_n_q_x_b[i] = q_x_b[i]*t_dot_n;
      t_dot_n_p_x_b[i] = p_x_b[i]*t_dot_n;
      n_q_x_b_dot_t[i] = q_x_b_dot_t*i_n[i];
      n_p_x_b_dot_t[i] = p_x_b_dot_t*i_n[i];
      // Calculating vectors.
      o_vec_int [0][i] = (t_dot_n_t_x_b[i] + t_t_x_b_dot_n[i])*i_one_m_nu        + b_x_t[i]*t_dot_n + t_b_x_t_dot_n[i]; // checked
      o_vec_int [1][i] = (t_dot_n_q_x_b[i] + t_q_x_b_dot_n[i])*i_one_m_nu        - n_q_x_b_dot_t[i] + i_q[i]*b_x_t_dot_n; // checked
      o_vec_int [2][i] = (t_dot_n_p_x_b[i] + t_p_x_b_dot_n[i])*i_one_m_nu        - n_p_x_b_dot_t[i] + i_p[i]*b_x_t_dot_n; // checked
      o_vec_int [3][i] = (t_dot_n_t_x_b[i] + t_t_x_b_dot_n[i])*one_m_nu_1p5_a_sq; // checked
      o_vec_int [4][i] = (t_dot_n_q_x_b[i] + t_q_x_b_dot_n[i])*one_m_nu_1p5_a_sq - n_q_x_b_dot_t[i]*a_sq_3; // checked
      o_vec_int [5][i] = (t_dot_n_p_x_b[i] + t_p_x_b_dot_n[i])*one_m_nu_1p5_a_sq - n_p_x_b_dot_t[i]*a_sq_3; // checked
      o_vec_int [6][i] = - t_dot_n_q_x_b_dot_t_3*i_q[i]; // checked
      o_vec_int [7][i] = - t_dot_n_p_x_b_dot_t_3*i_p[i]; // checked
      o_vec_int [8][i] = - t_dot_n_q_x_b_dot_t_3*i_t[i]; // checked
      o_vec_int [9][i] = - t_dot_n_p_x_b_dot_t_3*i_t[i]; // checked
      o_vec_int[10][i] = - t_dot_n_p_x_b_dot_t_3*i_q[i] - t_dot_n_q_x_b_dot_t_3*i_p[i]; // checked
  }
}

double *vertex_force_linear_rectangle(double *i_sch, double i_vec_int[][3],
                                      double i_r, double i_s, double i_rs,
                                      double i_factor,
                                      double *o_force){
  // Calculates the force on a vertex.
  double f[11];
  int i, j;
  f [0] = i_sch [3] - i_s*i_sch [6] - i_r*i_sch [8] + i_rs*i_sch [0];
  f [1] = i_sch [5] - i_s*i_sch [7] - i_r*i_sch[10] + i_rs*i_sch [2];
  f [2] = i_sch [4] - i_s*i_sch [9] - i_r*i_sch [7] + i_rs*i_sch [1];
  f [3] = i_sch[19] - i_s*i_sch[14] - i_r*i_sch[15] + i_rs*i_sch[11];
  f [4] = i_sch[21] - i_s*i_sch[16] - i_r*i_sch[18] + i_rs*i_sch[13];
  f [5] = i_sch[20] - i_s*i_sch[17] - i_r*i_sch[16] + i_rs*i_sch[12];
  f [6] = i_sch[35] - i_s*i_sch[28] - i_r*i_sch[32] + i_rs*i_sch[22];
  f [7] = i_sch[36] - i_s*i_sch[31] - i_r*i_sch[27] + i_rs*i_sch[23];
  f [8] = i_sch[33] - i_s*i_sch[26] - i_r*i_sch[30] + i_rs*i_sch[25];
  f [9] = i_sch[37] - i_s*i_sch[29] - i_r*i_sch[26] + i_rs*i_sch[24];
  f[10] = i_sch[34] - i_s*i_sch[27] - i_r*i_sch[28] + i_rs*i_sch[19];
  for (i = 0; i < 11; i++){
    for (j = 0; j < 3; j++){
      o_force[j] += i_vec_int[i][j] * f[i];
    }
  }
  o_force[0] = o_force[0]*i_factor;
  o_force[1] = o_force[1]*i_factor;
  o_force[2] = o_force[2]*i_factor;
  return o_force;
}

double *integrals_linear_rectangle(double *i_r, double *i_s, double *i_y,
                                   double *i_p, double *i_q, double *i_t,
                                   double i_p_dot_t, double i_q_dot_t, double i_a_sq,
                                   double *o_sch){//, int i_num_integrals = 38){
  // Sign vector for quick evaluation of integrals via the dot product.
  static double signv[8] = {1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0};
  static const double third = 1.0/3.0;
  static const double two_third = 2.0/3.0;
  // Integrals.
  double // Single integrals.
         a0m1[8], b0m1[8], c0m1[8], a1m1[8], b1m1[8], c1m1[8],
         a01 [8], b01 [8], c01 [8], a11 [8], b11 [8],
         // Double integrals.
         d00m3[8], e00m3[8], f00m3[8], d01m3[8], d10m3[8], d11m3[8], e01m3[8], e10m3[8],
         e11m3[8], f01m3[8], f10m3[8], f11m3[8], d00m1[8], e00m1[8], f00m1[8], d01m1[8],
         e01m1[8], f01m1[8], d10m1[8], e10m1[8], f10m1[8], d11m1[8], e11m1[8], f11m1[8],
         d02m3[8], d20m3[8], e20m3[8], f20m3[8], d001 [8], e001 [8], f001 [8], d02m1[8],
         d20m1[8], e20m1[8], f20m1[8], d12m3[8], d21m3[8], e21m3[8], f21m3[8], d22m3[8],
         // Triple integrals
         h[38][8] // Windows is dumb -.- .
		 //h [i_num_integrals][8]
		 , h001m1[8], h010m1[8], h100m1[8], h021m3[8], h201m3[8], h000m1[8];
 /*
   y^2, r^2, s^2 coordinates
   y*(p dot t), y*(q dot t), r*(p dot t), r*(q dot t)
    r dot p   ,  r dot q   ,  r dot t
   (r dot p)^2, (r dot q)^2, (r dot t)^2
 */
 double y_sq[8], r_sq[8], s_sq[8],
        y_p_dot_t[8], y_q_dot_t[8], r_p_dot_t[8], s_q_dot_t[8],
        r_dot_p[8], r_dot_q[8], r_dot_t[8],
        r_dot_p_sq[8], r_dot_q_sq[8], r_dot_t_sq[8];
  /*
   ra = sqrt((r_vec dot r_vec) + a**2) from non-singular dislocation theory, there are 8 r_vec, thus 8 ra's.
   ra_sq = ra^2 (element-wise squaring)
   ra_c_o_3 = 1/3 * ra^3 (element-wise cubing and division)
  */
  double ra[8], ra_sq[8], ra_c_o_3[8];
  double r_vec[8][3];
  /*
   Auxiliary constants to reduce computational cost.
   2*(p dot t)
   2*(q dot t)
   1-(p dot t)
   1-(q dot t)
   1-(p dot t)^2
   1-(q dot t)^2
   1/(1-(p dot t)^2)
   1/(1-(q dot t)^2)
   1-(p dot t)^2-(q dot t)^2
   1/(1-(p dot t)^2-(q dot t)^2)
   1/3 * 1/(1-(p dot t)^2-(q dot t)^2)
   (p dot t)*(q dot t)
  */
  double                         two_p_dot_t,
                                 two_q_dot_t,
                               one_m_p_dot_t,
                               one_m_q_dot_t,
                               one_m_p_dot_t_sq,
                               one_m_q_dot_t_sq,
                         one_o_one_m_p_dot_t_sq,
                         one_o_one_m_q_dot_t_sq,
                  one_m_p_dot_t_sq_m_q_dot_t_sq,
            one_o_one_m_p_dot_t_sq_m_q_dot_t_sq,
        trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq,
        p_dot_t_q_dot_t;
  int i, j;
  // Auxiliary scalars.
  two_p_dot_t = 2.0*i_p_dot_t;
  two_q_dot_t = 2.0*i_q_dot_t;
  one_m_p_dot_t = 1.0 - i_p_dot_t;
  one_m_q_dot_t = 1.0 - i_q_dot_t;
  one_m_p_dot_t_sq = 1.0 - i_p_dot_t * i_p_dot_t;
  one_m_q_dot_t_sq = 1.0 - i_q_dot_t * i_q_dot_t;
  p_dot_t_q_dot_t  =       i_q_dot_t * i_p_dot_t;
  one_o_one_m_p_dot_t_sq = 1.0 / one_m_p_dot_t_sq;
  one_o_one_m_q_dot_t_sq = 1.0 / one_m_q_dot_t_sq;
  one_m_p_dot_t_sq_m_q_dot_t_sq = one_m_p_dot_t_sq + one_m_q_dot_t_sq - 1.0;
  one_o_one_m_p_dot_t_sq_m_q_dot_t_sq = 1.0/one_m_p_dot_t_sq_m_q_dot_t_sq;
  trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq = third * one_o_one_m_p_dot_t_sq_m_q_dot_t_sq;
  // Auxiliary arrays.
  // Calculate r for all 8 combinations (from x1 to x3, x4, x5, x6 and x2 to x3, x4, x5, x6).
  for (i = 0; i < 8; i++){
    for (j = 0; j < 3; j++){
      r_vec[i][j] = i_r[i] * i_p[j] + i_s[i] * i_q[j] + i_y[i] * i_t[j];
    }
  }
  for (i = 0; i < 8; i++){
    r_sq      [i] = i_r[i] * i_r[i];
    s_sq      [i] = i_s[i] * i_s[i];
    y_sq      [i] = i_y[i] * i_y[i];
    ra        [i] = sqrt(r_sq[i]+s_sq[i]+y_sq[i] + two_p_dot_t*i_r[i]*i_y[i] + two_q_dot_t*i_s[i]*i_y[i] + i_a_sq);
    ra_sq     [i] = ra   [i] * ra[i];
    ra_c_o_3  [i] = ra_sq[i] * ra[i]/3.0;
    r_dot_p   [i] = dot_product(r_vec[i], i_p, 3);
    r_dot_q   [i] = dot_product(r_vec[i], i_q, 3);
    r_dot_t   [i] = dot_product(r_vec[i], i_t, 3);
    r_dot_p_sq[i] = r_dot_p[i] * r_dot_p[i];
    r_dot_q_sq[i] = r_dot_q[i] * r_dot_q[i];
    r_dot_t_sq[i] = r_dot_t[i] * r_dot_t[i];
    y_p_dot_t [i] = i_y[i] * i_p_dot_t;
    y_q_dot_t [i] = i_y[i] * i_q_dot_t;
    r_p_dot_t [i] = i_r[i] * i_p_dot_t;
    s_q_dot_t [i] = i_s[i] * i_q_dot_t;
  }
  // Calculate integrals.
  for (i = 0; i < 8; i++){
    // Linear seed integrals.
    a0m1[i] = seed_single_integral(ra[i], r_dot_p[i]); // checked
    b0m1[i] = seed_single_integral(ra[i], r_dot_q[i]); // checked
    c0m1[i] = seed_single_integral(ra[i], r_dot_t[i]); // checked
    // Linear integrals.
    // result = 0.5*(a*b + (c-d)*e)
    a01 [i] = integral_type_1(ra[i], r_dot_p[i], ra_sq[i], r_dot_p_sq[i], a0m1[i]); // checked
    b01 [i] = integral_type_1(ra[i], r_dot_q[i], ra_sq[i], r_dot_q_sq[i], b0m1[i]); // checked
    c01 [i] = integral_type_1(ra[i], r_dot_t[i], ra_sq[i], r_dot_t_sq[i], c0m1[i]); // checked
    // type_2 = a + b*c
    a1m1[i] = integral_type_2(ra[i]      , y_p_dot_t[i]             , -a0m1[i]); // checked
    b1m1[i] = integral_type_2(ra[i]      , y_q_dot_t[i]             , -b0m1[i]); // checked
    c1m1[i] = integral_type_2(ra[i]      , r_p_dot_t[i]+s_q_dot_t[i], -c0m1[i]); // checked
    a11 [i] = integral_type_2(ra_c_o_3[i], y_p_dot_t[i]             , -a01 [i]); // checked
    b11 [i] = integral_type_2(ra_c_o_3[i], y_q_dot_t[i]             , -b01 [i]); // checked
    // Double seed integrals.
    d00m3[i] = seed_double_integral(1.0, ra[i], r_dot_p[i], r_dot_q[i], 0.0,
                                    1.0             , i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, y_sq[i]); // checked
    e00m3[i] = seed_double_integral(one_m_p_dot_t, ra[i], i_r[i], i_y[i], s_q_dot_t[i],
                                    one_m_p_dot_t_sq, i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, s_sq[i]); // checked
    f00m3[i] = seed_double_integral(one_m_q_dot_t, ra[i], i_s[i], i_y[i], r_p_dot_t[i],
                                    one_m_q_dot_t_sq, i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, r_sq[i]); // checked
    // Double integrals.
    // type_2 = a + b*c
    d01m3[i] = integral_type_2(-a0m1[i], -y_q_dot_t[i], d00m3[i]); // checked
    d10m3[i] = integral_type_2(-b0m1[i], -y_p_dot_t[i], d00m3[i]); // checked
    // type_3 = a + b*c + d*e
    e01m3[i] = -one_o_one_m_p_dot_t_sq * integral_type_3(a0m1[i], s_q_dot_t[i], e00m3[i], -i_p_dot_t, c0m1[i]); // checked
    f01m3[i] = -one_o_one_m_q_dot_t_sq * integral_type_3(b0m1[i], r_p_dot_t[i], f00m3[i], -i_q_dot_t, c0m1[i]); // checked
    // type_2 = a + b*c
    e10m3[i] = integral_type_2(-c0m1[i], -i_p_dot_t, e01m3[i]); // checked
    f10m3[i] = integral_type_2(-c0m1[i], -i_q_dot_t, f01m3[i]); // checked
    // type_6 = a*b + c*d + (e+f)*g + h*i + j*k
    d00m1[i] = integral_type_6(i_r[i], b0m1[i], i_s[i], a0m1[i], i_a_sq, y_sq[i], -d00m3[i], -y_p_dot_t[i], d10m3[i], -y_q_dot_t[i], d01m3[i]); // checked
    // type_7 = a*b + c*d + (e+f)*g + h*i
    e00m1[i] = integral_type_7(i_r[i], c0m1[i], i_y[i], a0m1[i], i_a_sq, s_sq[i], -e00m3[i], -s_q_dot_t[i], e01m3[i]); // checked
    f00m1[i] = integral_type_7(i_s[i], c0m1[i], i_y[i], b0m1[i], i_a_sq, r_sq[i], -f00m3[i], -r_p_dot_t[i], f01m3[i]); // checked
    // type_2 = a + b*c
    d11m3[i] = integral_type_2(-a1m1[i], -y_q_dot_t[i], d10m3[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    e11m3[i] = one_o_one_m_p_dot_t_sq * integral_type_4(-a1m1[i], r_p_dot_t[i], c0m1[i], -i_p_dot_t, e00m1[i], -s_q_dot_t[i], e10m3[i]); // checked
    f11m3[i] = one_o_one_m_q_dot_t_sq * integral_type_4(-b1m1[i], s_q_dot_t[i], c0m1[i], -i_q_dot_t, f00m1[i], -r_p_dot_t[i], f10m3[i]); // checked
    // type_3 = a + b*c + d*e
    d20m3[i] = integral_type_3(d00m1[i], -y_p_dot_t[i], d10m3[i], -i_r[i], b0m1[i]); // checked
    d02m3[i] = integral_type_3(d00m1[i], -y_q_dot_t[i], d01m3[i], -i_s[i], a0m1[i]); // checked
    e20m3[i] = integral_type_3(e00m1[i], -i_p_dot_t     , e11m3[i], -i_r[i], c0m1[i]); // checked
    f20m3[i] = integral_type_3(f00m1[i], -i_q_dot_t     , f11m3[i], -i_s[i], c0m1[i]); // checked
    // type_2 = a + b*c
    d01m1[i] = integral_type_2(a01[i], -y_q_dot_t[i], d00m1[i]); // checked
    d10m1[i] = integral_type_2(b01[i], -y_p_dot_t[i], d00m1[i]); // checked
    e01m1[i] = one_o_one_m_p_dot_t_sq * integral_type_3(a01[i], -s_q_dot_t[i], e00m1[i], -i_p_dot_t, c01[i]); // checked
    f01m1[i] = one_o_one_m_q_dot_t_sq * integral_type_3(b01[i], -r_p_dot_t[i], f00m1[i], -i_q_dot_t, c01[i]); // checked
    // type_2 = a + b*c
    e10m1[i] = integral_type_2(c01[i], -i_p_dot_t, e01m1[i]); // checked
    f10m1[i] = integral_type_2(c01[i], -i_q_dot_t, f01m1[i]); // checked
    // type_6 = a*b + c*d + (e+f)*g + h*i + j*k
    d001[i]  = third * integral_type_6(i_r[i], b01[i], i_s[i], a01[i], y_sq[i], i_a_sq, d00m1[i], y_p_dot_t[i], d10m1[i], y_q_dot_t[i], d01m1[i]); // checked
    // type_7 = a*b + c*d + (e+f)*g + h*i
    e001[i]  = third * integral_type_7(i_r[i], c01[i], i_y[i], a01[i], s_sq[i], i_a_sq, e00m1[i], s_q_dot_t[i], e01m1[i]); // checked
    f001[i]  = third * integral_type_7(i_s[i], c01[i], i_y[i], b01[i], r_sq[i], i_a_sq, f00m1[i], r_p_dot_t[i], f01m1[i]); // checked
    // type_2 = a + b*c
    d11m1[i] = integral_type_2(a11[i], -y_q_dot_t[i], d10m1[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    e11m1[i] = one_o_one_m_p_dot_t_sq * integral_type_4(a11[i], -r_p_dot_t[i], c01[i], i_p_dot_t , e001[i], -s_q_dot_t[i], e10m1[i]); // checked
    f11m1[i] = one_o_one_m_q_dot_t_sq * integral_type_4(b11[i], -s_q_dot_t[i], c01[i], i_q_dot_t , f001[i], -r_p_dot_t[i], f10m1[i]); // checked
    // type_3 = a + b*c + d*e
    d02m1[i] = integral_type_3(-d001[i], -y_q_dot_t[i], d01m1[i],  i_s[i], a01 [i]); // checked
    d20m1[i] = integral_type_3(-d001[i], -y_p_dot_t[i], d10m1[i],  i_r[i], b01 [i]); // checked
    e20m1[i] = integral_type_3(-e001[i], -i_p_dot_t     , e11m1[i],  i_r[i], c01 [i]); // checked
    f20m1[i] = integral_type_3(-f001[i], -i_q_dot_t     , f11m1[i],  i_s[i], c01 [i]); // checked
    d12m3[i] = integral_type_3(d10m1[i], -y_q_dot_t[i], d11m3[i], -i_s[i], a1m1[i]); // checked
    d21m3[i] = integral_type_3(d01m1[i], -y_p_dot_t[i], d11m3[i], -i_r[i], b1m1[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    e21m3[i] = one_o_one_m_p_dot_t_sq * integral_type_5(e01m1[i], y_p_dot_t[i], a1m1[i], -i_p_dot_t, e10m1[i], p_dot_t_q_dot_t*i_s[i], e11m3[i], -i_r[i], c1m1[i]); // checked
    f21m3[i] = one_o_one_m_q_dot_t_sq * integral_type_5(f01m1[i], y_q_dot_t[i], b1m1[i], -i_q_dot_t, f10m1[i], p_dot_t_q_dot_t*i_r[i], f11m3[i], -i_s[i], c1m1[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    // type_8 = a*b + c*d + e*f + g*h
    d22m3[i] = integral_type_5(-d001[i], p_dot_t_q_dot_t*y_sq[i], d11m3[i], i_r[i], b01[i], -i_r[i]*i_s[i], ra[i], i_s[i], a01[i])
             + i_y[i] * integral_type_8(-i_p_dot_t, d10m1[i], -i_q_dot_t, d01m1[i], i_p_dot_t*i_s[i], a1m1[i], i_q_dot_t*i_r[i], b1m1[i]); // checked
    // Triple integrals
    // type_3 = a + b*c + d*e
    h[0][i] = -one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_3(d00m1[i], -i_q_dot_t, e00m1[i], -i_p_dot_t, f00m1[i]); // checked
    // type_2 = a + b*c
    h[1][i] = integral_type_2(-f00m1[i], -i_p_dot_t, h[0][i]); // checked
    h[2][i] = integral_type_2(-e00m1[i], -i_q_dot_t, h[0][i]); // checked
    // type_3 = a + b*c + d*e
    h001m1[i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_3(d001[i], -i_q_dot_t, e001[i], -i_p_dot_t, f001[i]); // checked
    // type_2 = a + b*c
    h100m1[i] = integral_type_2(f001[i], -i_p_dot_t, h001m1[i]); // checked
    h010m1[i] = integral_type_2(e001[i], -i_q_dot_t, h001m1[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[3][i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(-d11m1[i], r_p_dot_t[i], f10m1[i], -i_p_dot_t, h010m1[i], s_q_dot_t[i], e10m1[i], -i_q_dot_t, h100m1[i]); // checked
    // type_3 = a + b*c + d*e
    h[4][i] = integral_type_3(h010m1[i], -i_p_dot_t, h[3][i], -i_r[i], f10m1[i]); // checked
    h[5][i] = integral_type_3(h100m1[i], -i_q_dot_t, h[3][i], -i_s[i], e10m1[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h021m3[i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d02m1[i], -two_q_dot_t, h010m1[i], i_p_dot_t, f20m1[i], i_q_dot_t*s_sq[i], e00m1[i]); // checked
    h201m3[i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d20m1[i], -two_p_dot_t, h100m1[i], i_q_dot_t, e20m1[i], i_p_dot_t*r_sq[i], f00m1[i]);
    // type_8 = a*b + c*d + e*f + g*h
    h000m1[i] = 0.5 * integral_type_8(-i_a_sq, 0.0, i_y[i], d00m1[i], i_s[i], e00m1[i], i_r[i], f00m1[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h[6][i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d10m1[i], r_p_dot_t[i], f00m1[i], i_q_dot_t, e10m1[i], -i_p_dot_t, h000m1[i]); // checked
    h[8][i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d01m1[i], s_q_dot_t[i], e00m1[i], i_p_dot_t, f10m1[i], -i_q_dot_t, h000m1[i]); // checked
    // type_2 = a + b*c
    h[7][i] = integral_type_2(-e10m1[i], -i_q_dot_t, h[6][i]); // checked
    // type_3 = a + b*c + d*e
    h[9] [i] = integral_type_3(h000m1[i], -i_p_dot_t, h[6][i], -i_r[i], f00m1[i]); // checked
    h[10][i] = integral_type_3(h000m1[i], -i_q_dot_t, h[8][i], -i_s[i], e00m1[i]); // checked
    h[11][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_3(-d00m3[i], i_p_dot_t, f00m3[i], i_q_dot_t, e00m3[i]); // checked
    // type_2 = a + b*c
    h[12][i] = integral_type_2(-third * f00m3[i], -i_p_dot_t, h[11][i]); // checked
    h[13][i] = integral_type_2(-third * e00m3[i], -i_q_dot_t, h[11][i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h[14][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d10m3[i], r_p_dot_t[i], f00m3[i], i_q_dot_t, e10m3[i], -i_p_dot_t, 0.0); // checked
    h[15][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d01m3[i], s_q_dot_t[i], e00m3[i], i_p_dot_t, f10m3[i], -i_q_dot_t, 0.0); // checked
    // type_2 = a + b*c
    h[16][i] = integral_type_2(-third * f10m3[i], -i_p_dot_t, h[15][i]); // checked
    // type_3 = a + b*c + d*e
    h[17][i] = integral_type_3(third * 0.0, -third*i_r[i], f00m3[i], -i_p_dot_t, h[14][i]); // checked
    h[18][i] = integral_type_3(third * 0.0, -third*i_s[i], e00m3[i], -i_q_dot_t, h[15][i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[19][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(-d11m3[i], r_p_dot_t[i], f10m3[i], -i_p_dot_t, h[2][i], s_q_dot_t[i], e10m3[i], -i_q_dot_t, h[1][i]); // checked
    // type_3 = a + b*c + d*e
    h[20][i] = integral_type_3(third * h[2][i], -third*i_r[i], f10m3[i], -i_p_dot_t, h[19][i]); // checked
    h[21][i] = integral_type_3(third * h[1][i], -third*i_s[i], e10m3[i], -i_q_dot_t, h[19][i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h[23][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d20m3[i], -two_p_dot_t, h[1][i], i_q_dot_t, e20m3[i], i_p_dot_t*r_sq[i], f00m3[i]); // checked
    h[22][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d02m3[i], -two_q_dot_t, h[2][i], i_p_dot_t, f20m3[i], i_q_dot_t*s_sq[i], e00m3[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[25][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(h[2][i], s_q_dot_t[i], e01m3[i], -i_y[i], d01m3[i], i_p_dot_t, f11m3[i], -i_q_dot_t, h[0][i]); // checked
    h[24][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(h[1][i], r_p_dot_t[i], f01m3[i], -i_y[i], d10m3[i], i_q_dot_t, e11m3[i], -i_p_dot_t, h[0][i]); // checked
    // type_9 = a + b*c + d*e + f*g + h*i + j*k
    h[26][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_9(h[7][i], r_p_dot_t[i], f11m3[i], -i_y[i], d11m3[i], s_q_dot_t[i], e11m3[i], -i_p_dot_t, h[8][i], -i_q_dot_t, h[6][i]); // checked
    // type_3 = a + b*c + d*e
    h[27][i] = integral_type_3(third * h[8][i], -i_p_dot_t, h[26][i], -third*i_r[i], f11m3[i]); // checked
    h[28][i] = integral_type_3(third * h[6][i], -i_q_dot_t, h[26][i], -third*i_s[i], e11m3[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[29][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(h[9] [i], -two_p_dot_t, h[6][i], i_q_dot_t, e21m3[i], i_p_dot_t*r_sq[i], f01m3[i], -i_y[i], d20m3[i]); // checked
    h[30][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(h[10][i], -two_q_dot_t, h[8][i], i_p_dot_t, f21m3[i], i_q_dot_t*s_sq[i], e01m3[i], -i_y[i], d02m3[i]); // checked
    // type_3 = a + b*c + d*e
    h[31][i] = integral_type_3(two_third * h[6][i], -i_p_dot_t, h[29][i], -third*r_sq[i], f01m3[i]); // checked
    h[32][i] = integral_type_3(two_third * h[8][i], -i_q_dot_t, h[30][i], -third*s_sq[i], e01m3[i]); // checked
    // type_10 = a + b*c + d*e + f*g + h*i + j*k
    h[33][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_10(h[5][i], -two_q_dot_t, h[3][i], -i_y[i], d12m3[i], r_p_dot_t[i], f21m3[i], -i_p_dot_t, h021m3[i], i_q_dot_t*s_sq[i], e11m3[i]); // checked
    h[37][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_10(h[4][i], -two_p_dot_t, h[3][i], -i_y[i], d21m3[i], s_q_dot_t[i], e21m3[i], -i_q_dot_t, h201m3[i], i_p_dot_t*r_sq[i], f11m3[i]); // checked
    // type_11 = a + b*c + d*e + f*g + h*i
    h[34][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_11(-d22m3[i], i_p_dot_t*r_sq[i], f20m3[i], i_q_dot_t*s_sq[i], e20m3[i], -i_p_dot_t*2.0, h[5][i], -i_q_dot_t*2.0, h[4][i]);
    // type_12 = a*b + c*d + e*f + g*h + i*j + k*l*m;
    h[36][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_12(2.0*one_m_q_dot_t_sq, h[3][i], -i_p_dot_t, h[4][i], p_dot_t_q_dot_t, h201m3[i], y_p_dot_t[i], d21m3[i], -p_dot_t_q_dot_t*i_s[i], e21m3[i], -one_m_q_dot_t_sq, r_sq[i], f11m3[i]); // checked
    h[35][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_12(2.0*one_m_p_dot_t_sq, h[3][i], -i_q_dot_t, h[5][i], p_dot_t_q_dot_t, h021m3[i], y_q_dot_t[i], d12m3[i], -p_dot_t_q_dot_t*i_r[i], f21m3[i], -one_m_p_dot_t_sq, s_sq[i], e11m3[i]); // checked
  }
  // Evaluating the integrals.
  for (i = 0; i<38; i++){
    o_sch[i] = dot_product(h[i], signv, 8);
  }
  return o_sch;
}

void compute_forces_linear_rectangle(double *i_sch, double i_vec_int[][3],
                                     double *i_rp, double *i_sp, double i_factor,
                                     double *o_nodal_force[4], double *o_total_force){
  int i;
  // Calculating nodal forces
  // x3.
  vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[1], i_sp[1], i_rp[1]*i_sp[1],  i_factor, o_nodal_force[0]);
  // x4.
  vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[0], i_sp[1], i_rp[0]*i_sp[1], -i_factor, o_nodal_force[1]);
  // x5.
  vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[1], i_sp[0], i_rp[1]*i_sp[0], -i_factor, o_nodal_force[2]);
  // x5
  vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[0], i_sp[0], i_rp[0]*i_sp[0],  i_factor, o_nodal_force[3]);
  for (i = 0; i < 4; i++){
    o_total_force[0] += o_nodal_force[i][0];
    o_total_force[1] += o_nodal_force[i][1];
    o_total_force[2] += o_nodal_force[i][2];
  }
}

void init_force(double *nodal_force[n_nodes], double *total_force){
  // Sets forces to zero.
  int i, j;
  for (i = 0; i < n_nodes; i++){
    for (j = 0; j < 3; j++){
      nodal_force[i][j] = 0.0;
    }
  }
  for (i = 0; i < 3; i++){
    total_force[i] = 0.0;
  }
}

void add_force(double *p_nodal_force[4], double *p_total_force, double *nodal_force[4], double *total_force){
  // Adds forces for averaging purposes later on.
  int i, j;
  for (i = 0; i < 4; i++){
    for (j = 0; j < 3; j++){
      nodal_force[i][j] += p_nodal_force[i][j];
    }
  }
  for (i = 0; i < 3; i++){
    total_force[i] += p_total_force[i];
  }
}

void mean_force(double *nodal_force[4], double *total_force, int n_samples){
  // Sets forces to zero.
  int i, j;
  for (i = 0; i < 4; i++){
    for (j = 0; j < 3; j++){
      nodal_force[i][j] = nodal_force[i][j]/n_samples;
    }
  }
  for (i = 0; i < 3; i++){
    total_force[i] = total_force[i]/n_samples;
  }
}

void nodal_surface_force_linear_rectangle(double *x1, double *x2, double *x3, double *x4, double *x5, double *x6, double *b, double *p, double *q, double *n, double p_norm, double q_norm, double mu, double nu, double a, double a_sq, double one_m_nu, double factor, double *nodal_force[4], double *total_force){
  // Vertices of the dislocation segment and a linear rectangular surface.
  /*
          x5 -------------------- x6
             |                  |                 t
             |                  |               ---->           b
             |                  |           x1 -------- x2    ^>
             |                  |             /        |     /
      ^      |                  |            /         \
    q |      |                  |          ...        ...
      |   X3 -------------------- X4

          ---->
            p
  */
  // Characteristic vectors.
  double t[3];
  // Basis vectors (unitary).
  double p_x_t[3], q_x_t[3];
  // Limits of the distance vector from the plane to dislocation line segment.
  // r_lim[0][] = vector from x1 to x3, r_lim[1][] = vector from x2 to x6.
  double r_lim[2][3];
  // r, s limits
  double rp[2], sp[2];
  //  y, r, s coordinates
  double y[8], r[8], s[8];
  // Vectors for the integrals.
  double vec_int[11][3];
  // Scalar value of the integrals.
  double sch[38];
  //  Auxiliary constants to reduce computational cost.
  //  p dot t, q dot t
  double p_dot_t, q_dot_t;
  int i;
  // Set forces to zero.
  init_force(nodal_force, total_force);
  // Build unit vectors t.
  init_vector (x1, x2, 3, t);
  // Dot products.
  p_dot_t = dot_product(p, t, 3);
  q_dot_t = dot_product(q, t, 3);
  //*******************WARNING*******************//
  // This formulation assumes x3-x6 is diagonal! //
  //*******************WARNING*******************//
  cross_product(p, t, p_x_t);
  cross_product(q, t, q_x_t);
  // Vectors between x3 and x1, and x6 and x2.
  for (i = 0; i < 3; i++){
    r_lim[0][i] = x3[i] - x1[i];
    r_lim[1][i] = x6[i] - x2[i];
  }
  // Integral bounds for y, r, s.
  for (i = 0; i < 2; i++){
    rp[i] = init_point(r_lim[i], q_x_t, p, q_x_t, 3);
    sp[i] = init_point(r_lim[i], p_x_t, q, p_x_t, 3);
  }
  // Assign coordinates for the evaluation of the integrals.
  y[0] = y[2] = y[4] = y[6] = init_point(r_lim[1], n    , t, n    , 3);
  y[1] = y[3] = y[5] = y[7] = init_point(r_lim[0], n    , t, n    , 3);
  r[0] = r[1] = r[2] = r[3] = rp[1];
  r[4] = r[5] = r[6] = r[7] = rp[0];
  s[0] = s[1] = s[4] = s[5] = sp[1];
  s[2] = s[3] = s[6] = s[7] = sp[0];
  // Calculate vectors for integrals.
  integral_vector(p, q, b, t, n, one_m_nu, a_sq, vec_int);
  // Calculate integrals.
  integrals_linear_rectangle(r, s, y, p, q, t, p_dot_t, q_dot_t, a_sq, sch);//, 38);
  // Calculate nodal forces.
  compute_forces_linear_rectangle(sch, vec_int, rp, sp, factor, nodal_force, total_force);
  //printf("total_force[x, y, z] = [%2.14f, %2.14f, %2.14f]\n", total_force[0], total_force[1], total_force[2]);
}

void nodal_surface_force_linear_rectangle_special(double *x1, double *x2, double *x3, double *x4, double *x5, double *x6, double *b, double *t, double *p, double *q, double *n, double p_norm, double q_norm, double mu, double nu, double a, double a_sq, double one_m_nu, double factor, double *nodal_force[n_nodes], double *total_force){
  /*
    Forces
    nodal_force[0][] = F_x3[x, y, z], nodal_force[1][] = F_x4[x, y, z],
    nodal_force[2][] = F_x5[x, y, z], nodal_force[3][] = F_x6[x, y, z]
    total_force[x, y, z] = F_x3[x, y, z] + F_x4[x, y, z] + F_x5[x, y, z] + F_x6[x, y, z]
  */
  // Modulus of p and q.
  double rot_centre[3], rot_x1[3], rot_x2[3];
  double t_x_n[3], mag_t_x_n, p_total_force[3], *p_nodal_force[n_nodes];
  int rotation;
  rotation = 4;
  // Initialise force to zero.
  init_force(nodal_force, total_force);

  for (int i = 0; i < n_nodes; i++){
    p_nodal_force[i] = (double *) malloc(3 * sizeof(double));
  }
  double angle = pi/180.;
  cross_product(t, n, t_x_n);
  mag_t_x_n = sqrt(dot_product(t_x_n, t_x_n, 3));
  for (int i = 0; i < 3; i++){
    // Halfway between x1 and x2. x1 + (x2-x1)/2
    rot_centre[i] = 0.5*(x1[i] + x2[i]);
    t_x_n[i] = t_x_n[i]/mag_t_x_n;
  }
  for (int i = 1; i < rotation + 1; i++){
    arbitrary_rotation_matrix_3d(i*angle, rot_centre, t_x_n, x1, rot_x1);
    arbitrary_rotation_matrix_3d(i*angle, rot_centre, t_x_n, x2, rot_x2);
    nodal_surface_force_linear_rectangle(rot_x1, rot_x2, x3, x4, x5, x6, b, p, q, n, p_norm, q_norm, mu, nu, a, a_sq, one_m_nu, factor, p_nodal_force, p_total_force);
    add_force(p_nodal_force, p_total_force, nodal_force, total_force);
    arbitrary_rotation_matrix_3d(-i*angle, rot_centre, t_x_n, x1, rot_x1);
    arbitrary_rotation_matrix_3d(-i*angle, rot_centre, t_x_n, x2, rot_x2);
    nodal_surface_force_linear_rectangle(rot_x1, rot_x2, x3, x4, x5, x6, b, p, q, n, p_norm, q_norm, mu, nu, a, a_sq, one_m_nu, factor, p_nodal_force, p_total_force);
    add_force(p_nodal_force, p_total_force, nodal_force, total_force);
  }
  mean_force(nodal_force, total_force, rotation*2);
  for (int i = 0; i < n_nodes; i++){
    free(p_nodal_force[i]);
  }
}

void main_nodal_surface_force_linear_rectangle(double *x1, double *x2, double *x3, double *x4, double *x5, double *x6, double *b, double mu, double nu, double a, double a_sq, double one_m_nu, double factor, double *nodal_force[4], double *total_force){
  /*
    Forces
    nodal_force[0][] = F_x3[x, y, z], nodal_force[1][] = F_x4[x, y, z],
    nodal_force[2][] = F_x5[x, y, z], nodal_force[3][] = F_x6[x, y, z]
    total_force[x, y, z] = F_x3[x, y, z] + F_x4[x, y, z] + F_x5[x, y, z] + F_x6[x, y, z]
  */
  // Characteristic vectors.
  double p[3], q[3], n[3], t[3];
  // Modulus of p and q.
  double p_norm, q_norm, t_dot_n;
  double l_factor;
  double rot_centre[3], rot_x1[3], rot_x2[3];
  double t_x_n[3], mag_t_x_n, p_total_force[3], *p_nodal_force[4];
  double angle = 0.01*pi/180.;
  int rotation;
  int i, j;
  rotation = 3;
  // Vectors
  init_vector(x1, x2, 3, t);
  init_vector2(x3, x4, 3, p, &p_norm);
  init_vector2(x3, x5, 3, q, &q_norm);
  cross_product(p, q, n);
  normalise_vector(n, 3, n);
  t_dot_n = dot_product(t, n, 3);
  // Auxiliary variables.
  l_factor = factor/p_norm/q_norm;
  if(t_dot_n != 0.0){
    nodal_surface_force_linear_rectangle(x1, x2, x3, x4, x5, x6, b, p, q, n, p_norm, q_norm, mu, nu, a, a_sq, one_m_nu, l_factor, nodal_force, total_force);
  }
  else{
	// Initialise force to zero.
	init_force(nodal_force, total_force);
    for (i = 0; i < 4; i++){
      p_nodal_force[i] = (double *) malloc(3*sizeof(double));
    }
    cross_product(t, n, t_x_n);
    mag_t_x_n = sqrt(dot_product(t_x_n, t_x_n, 3));
    for (i = 0; i < 3; i++){
      // Halfway between x1 and x2. x1 + (x2-x1)/2
      rot_centre[i] = 0.5*(x1[i] + x2[i]);
      t_x_n[i] = t_x_n[i]/mag_t_x_n;
    }
    //FILE *fp;
    //fp = fopen("./tests/test2.txt", "w");
    for (i = 0; i < rotation; i++){
      j = i+1;
      arbitrary_rotation_matrix_3d(j*angle, rot_centre, t_x_n, x1, rot_x1);
      arbitrary_rotation_matrix_3d(j*angle, rot_centre, t_x_n, x2, rot_x2);
      nodal_surface_force_linear_rectangle(rot_x1, rot_x2, x3, x4, x5, x6, b, p, q, n, p_norm, q_norm, mu, nu, a, a_sq, one_m_nu, l_factor, p_nodal_force, p_total_force);
      //printf("theta  = %f\nrot_x1 = [%f, %f, %f]\nrot_x2 = [%f, %f, %f]\nx1 = [%f, %f, %f]\nx2 = [%f, %f, %f]\n", j*angle, rot_x1[0], rot_x1[1], rot_x1[2], rot_x2[0], rot_x2[1], rot_x2[2], x1[0], x1[1], x1[2], x2[0], x2[1], x2[2]);
      //fprintf(fp, "%3f %2.14f %2.14f %2.14f\n", j*angle*180./pi, p_total_force[0], p_total_force[1], p_total_force[2]);
      add_force(p_nodal_force, p_total_force, nodal_force, total_force);
      arbitrary_rotation_matrix_3d(-j*angle, rot_centre, t_x_n, x1, rot_x1);
      arbitrary_rotation_matrix_3d(-j*angle, rot_centre, t_x_n, x2, rot_x2);
      nodal_surface_force_linear_rectangle(rot_x1, rot_x2, x3, x4, x5, x6, b, p, q, n, p_norm, q_norm, mu, nu, a, a_sq, one_m_nu, l_factor, p_nodal_force, p_total_force);
      //printf("theta  = %f\nrot_x1 = [%f, %f, %f]\nrot_x2 = [%f, %f, %f]\nx1 = [%f, %f, %f]\nx2 = [%f, %f, %f]\n", j*angle, rot_x1[0], rot_x1[1], rot_x1[2], rot_x2[0], rot_x2[1], rot_x2[2], x1[0], x1[1], x1[2], x2[0], x2[1], x2[2]);
      //fprintf(fp, "%3.6f %2.14f %2.14f %2.14f\n", -j*angle*180./pi, p_total_force[0], p_total_force[1], p_total_force[2]);
      add_force(p_nodal_force, p_total_force, nodal_force, total_force);
    }
    //fclose(fp);
    mean_force(nodal_force, total_force, rotation*2);
    for (i = 0; i < 4; i++){
      free(p_nodal_force[i]);
    }
  }
  //printf("total_force[x, y, z] = [%2.14f, %2.14f, %2.14f]\n", total_force[0], total_force[1], total_force[2]);
}

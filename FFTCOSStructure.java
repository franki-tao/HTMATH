package 结构性存款;

import org.apache.commons.math3.complex.Complex;
import org.apache.spark.internal.config.R;
import 结构性存款.product.product;

import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

public class FFTCOSStructure {
    protected final SNUtil util = new SNUtil();
    protected final double pi = Math.PI;
    public double L; //积分区间长度
    public int[] obs_list; //敲出观察日
    public double q; //红利
    public int n; //傅里叶展开项数
    public int N_steps; //时间步长
    public double HU; //敲出价
    public double HL; //敲入价
    public double[] c; //积分区间参数
    public double S_0; //期初价格
    public double T; //期限
    public double r; //无风险利率
    public double d; //红利
    public double strike; //行权价
    public Map<String, Double> param_dict; //参数字典
    public double mu_f; //息票率
    public double P; //本金
    public String model; //可选模型
    public product pro = new product(); //结构性产品


    public FFTCOSStructure() {
    }

    protected double[] COS_DUO_1(int n, int N_step, double[] c, int cp, double S_0, double HU, double HL, double T, double r, double d,
                                 double strike, Map<String, Double> param_dict, int[] UO_list, double mu_f, double P, String model) {
        if (pro.ID == 2 || pro.ID == 3 || pro.ID == 11) {
            return new double[]{0, 0, 0};
        }
        Complex iutil = new Complex(0, 1);
        double dt = T / N_step;
        double N = Math.pow(2, n);
        double x = Math.log(S_0 / strike);
        double hu = Math.log(HU / strike);
        double hl = Math.log(HL / strike);
        int idx_temp = 0;
        if (pro.ID == 9) {
            mu_f = pro.kd_dict.get("mu_f");
        }
        double cash = (1 + mu_f * T) * P;
        if (pro.ID == 4 || pro.ID == 7 || pro.ID == 13) {
            cash = (1 + pro.kd_dict.get("mu_f") * T) * P;
        } else if (pro.ID == 10) {
            cash = (1 + pro.kd_dict.get("mu_f") * P * obs_list.length);
        } else if (pro.ID == 12) {
            obs_list = pro.obs_list;
            cash = 1 + pro.kd_dict.get("mu_f");
        } else if (pro.ID == 16) {
            cash = (1 + pro.kd_dict.get("min_rate") * T - pro.kd_dict.get("opt_val")) * P;
        } else if (pro.ID == 17) {
            cash = (1 + pro.kd_dict.get("min_rate") * T) * P;
        } else if (pro.ID == 18) {
            cash = (1 + pro.kd_dict.get("mu_f") * T - pro.kd_dict.get("opt_val")) * P;
        } else if (pro.ID == 19) {
            cash = (1 + pro.kd_dict.get("mu_f") * T) * P;
        }

        double a = c[0] + x - L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double b = c[0] + x + L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double[][] Grid_i = new double[(int) N][1];
        Complex[] aux = new Complex[(int) N];
        for (int i = 0; i < Grid_i.length; i++) {
            Grid_i[i][0] = i;
            aux[i] = new Complex(pi * i / (b - a));
        }
        double[][] V = util.Coefficient_V(Grid_i, new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b}, cash);

        if (pro.ID == 5 || pro.ID == 6 || pro.ID == 7 || pro.ID == 13 || pro.ID == 15) {
            idx_temp = pro.KO_arr.length - 1;
            HU = pro.KO_arr[idx_temp];
            V = util.Coefficient_V(Grid_i, new double[]{hl}, new double[]{Math.log(HU / strike)}, new double[]{a},
                    new double[]{b}, cash);
        } else if (pro.ID == 24) {
            idx_temp = pro.KL_arr.length - 1;
            HL = pro.KL_arr[idx_temp];
            V = util.Coefficient_V(Grid_i, new double[]{Math.log(HL / strike)}, new double[]{hu}, new double[]{a},
                    new double[]{b}, cash);
        }

        for (int i = N_step - 1; i > 0; i--) {
            if (this.contain(obs_list, i)) {
                if (pro.ID == 5 || pro.ID == 6 || pro.ID == 7 || pro.ID == 13 || pro.ID == 15) {
                    idx_temp -= 1;
                    HU = pro.KO_arr[idx_temp];
                    hu = Math.log(HU / strike);
                    V = util.C_value(new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b}, (int) N, V, dt, r,
                            d, q, model, param_dict);
                } else if (pro.ID == 24) {
                    HL = pro.KL_arr[i - 1];
                    hl = Math.log(HL / strike);
                    V = util.C_value(new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b}, (int) N, V, dt, r,
                            d, q, model, param_dict);
                } else {
                    V = util.C_value(new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                }
            } else {
                if (pro.ID == 23) {
                    V = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                } else if (pro.ID == 24) {
                    HL = pro.KL_arr[i - 1];
                    hl = Math.log(HL / strike);
                    V = util.C_value(new double[]{hl}, new double[]{b}, new double[]{a}, new double[]{b}, (int) N, V, dt, r,
                            d, q, model, param_dict);
                } else {
                    V = util.C_value(new double[]{hl}, new double[]{b}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                }
            }
        }


        Complex[] cf_value = util.process(aux, dt, r, d, q, param_dict, model);
        for (int i = 0; i < cf_value.length; i++) {
            cf_value[i] = cf_value[i].exp();
        }
        Complex[] pF = new Complex[cf_value.length];
        Complex[] dF = new Complex[cf_value.length];
        Complex[] gF = new Complex[cf_value.length];
        for (int i = 0; i < pF.length; i++) {
            pF[i] = cf_value[i].multiply(iutil.multiply(aux[i]).multiply(x - a).exp());
            dF[i] = pF[i].multiply(iutil).multiply(aux[i]);
            gF[i] = pF[i].multiply(iutil.multiply(-1).multiply(aux[i]).add(iutil.multiply(aux[i]).multiply(iutil).multiply(aux[i])));
        }
        pF[0] = pF[0].multiply(0.5);
        dF[0] = dF[0].multiply(0.5);
        gF[0] = gF[0].multiply(0.5);
        double price = 0d, delta = 0d, gamma = 0d;
        for (int i = 0; i < pF.length; i++) {
            price += pF[i].getReal() * V[i][0];
            delta += dF[i].getReal() * V[i][0] / S_0;
            gamma += gF[i].getReal() * V[i][0] / (S_0 * S_0);
        }
        price *= Math.exp(-r * dt);
        delta *= Math.exp(-r * dt);
        gamma *= Math.exp(-r * dt);
        return new double[]{price, delta, gamma};
    }

    protected double[] COS_UO_2(int n, int N_step, double[] c, int cp, double S_0, double HU, double T, double r, double d, double strike,
                                Map<String, Double> param_dict, int[] UO_list, double mu_f, double P, String model) {
        Complex iutil = new Complex(0, 1);
        double dt = T / N_step;
        double N = Math.pow(2, n);
        double x = Math.log(S_0 / strike);
        double hu = Math.log(HU / strike);
        double cash = (1 + mu_f * T) * P;
        int idx_temp = 0;
        int obs_list_len_temp = 0;
        if (pro.ID == 4 || pro.ID == 7 || pro.ID == 13) {
            idx_temp = pro.mu_O.length - 1;
            cash = (1 + pro.mu_O[idx_temp] * T) * P;
        } else if (pro.ID == 10 || pro.ID == 11) {
            obs_list_len_temp = obs_list.length;
            cash = (1 + pro.kd_dict.get("mu_f") * obs_list_len_temp) * P;
        } else if (pro.ID == 12) {
            obs_list = pro.obs_list;
            cash = (1 + pro.kd_dict.get("mu_f"));
        } else if (pro.ID == 16 || pro.ID == 18) {
            cash = (1 + pro.kd_dict.get("mu_O") * T - pro.kd_dict.get("opt_val")) * P;
        } else if (pro.ID == 17 || pro.ID == 19) {
            cash = (1+pro.kd_dict.get("mu_O")*T)*P;
        }
        double a = c[0] + x - L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double b = c[0] + x + L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double[][] Grid_i = new double[(int) N][1];
        Complex[] aux = new Complex[(int) N];
        for (int i = 0; i < Grid_i.length; i++) {
            Grid_i[i][0] = i;
            aux[i] = new Complex(pi * i / (b - a));
        }
        double[][] V = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, cash);
        if (pro.ID == 3) {
            double hu_2 = Math.log(pro.kd_dict.get("K_2") / strike);
            cash = (1 + pro.kd_dict.get("mu_O")) * P;
            V = util.Coefficient_V(Grid_i, new double[]{hu_2}, new double[]{b}, new double[]{a}, new double[]{b}, cash);
        } else if (pro.ID == 5 || pro.ID == 6 || pro.ID == 7 || pro.ID == 13 || pro.ID == 15) {
            idx_temp = pro.KO_arr.length - 1;
            HU = pro.KO_arr[idx_temp];
            hu = Math.log(HU / strike);
            V = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, cash);
        } else if (pro.ID == 8) {
            double[][] tp1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, cash);
            double[][] tp2 = util.Coefficient_V_kd8(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b});
            double[][] tp3 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, 1);
            for (int i = 0; i < tp1.length; i++) {
                for (int j = 0; j < tp1[0].length; j++) {
                    V[i][j] = tp1[i][j] + (tp2[i][j] - HU / S_0 * tp3[i][j]) * pro.kd_dict.get("prt_rate");
                }
            }
        } else if (pro.ID == 23) {
            HL = pro.kd_dict.get("HL");
            double hl = Math.log(HL / strike);
            V = util.Coefficient_V(Grid_i, new double[]{a}, new double[]{hl}, new double[]{a}, new double[]{b}, cash);
        }
        for (int m = N_step - 1; m > 0; m--) {
            if (this.contain(obs_list, m)) {
                if (pro.ID == 3) {
                    double Rb = (1 + pro.kd_dict.get("mu_O")) * P;
                    double hu_2 = Math.log(pro.kd_dict.get("K_2") / strike);
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu_2}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu_2}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else if (pro.ID == 4) {
                    idx_temp -= 1;
                    double Rb = (1 + pro.mu_O[idx_temp] * dt * m) * P;
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else if (pro.ID == 5 || pro.ID == 6) {
                    idx_temp -= 1;
                    double Rb = (1 + mu_f * dt * m) * P;
                    HU = pro.KO_arr[idx_temp];
                    hu = Math.log(HU / strike);
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else if (pro.ID == 7 || pro.ID == 13) {
                    idx_temp -= 1;
                    double Rb = (1 + pro.mu_O[idx_temp] * dt * m) * P;
                    HU = pro.KO_arr[idx_temp];
                    hu = Math.log(HU / strike);
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else if (pro.ID == 8) {
                    double Rb = (1 + mu_f * dt * m) * P;
                    double[][] tp1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] tp2 = util.Coefficient_V_kd8(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b});
                    double[][] tp3 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, 1);
                    double[][] G1 = new double[tp1.length][tp1[0].length];
                    for (int i = 0; i < G1.length; i++) {
                        for (int j = 0; j < G1[0].length; j++) {
                            G1[i][j] = tp1[i][j] + (tp2[i][j] - HU / S_0 * tp3[i][j]) * pro.kd_dict.get("prt_rate");
                        }
                    }
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);

                } else if (pro.ID == 10 || pro.ID == 11) {
                    obs_list_len_temp -= 1;
                    double Rb = (1 + pro.kd_dict.get("mu_f") * obs_list_len_temp) * P;
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else if (pro.ID == 12) {
                    double Rb = (1 + pro.kd_dict.get("mu_f")) * P;
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else if (pro.ID == 15) {
                    idx_temp -= 1;
                    double Rb = (1 + mu_f * dt * m) * P;
                    HU = pro.KO_arr[idx_temp];
                    hu = Math.log(HU / strike);
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else if (pro.ID == 16 || pro.ID == 18) {
                    double Rb = (1 + pro.kd_dict.get("mu_O") * dt * m - pro.kd_dict.get("opt_val")) * P;
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else if (pro.ID == 17 || pro.ID == 19) {
                    double Rb = (1 + pro.kd_dict.get("mu_O") * dt * m) * P;
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else if (pro.ID == 23) {
                    HL = pro.kd_dict.get("HL");
                    double hl = Math.log(HL / strike);
                    double Rb = (1 + mu_f * dt * m) * P;
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{a}, new double[]{hl}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{hl}, new double[]{b}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                } else {
                    double Rb = (1 + mu_f * dt * m) * P;
                    double[][] G1 = util.Coefficient_V(Grid_i, new double[]{hu}, new double[]{b}, new double[]{a}, new double[]{b}, Rb);
                    double[][] G2 = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                    V = this.add(G1, G2);
                }

            } else {
                V = util.C_value(new double[]{a}, new double[]{b}, new double[]{a}, new double[]{b},
                        (int) N, V, dt, r, d, q, model, param_dict);
            }
        }

        Complex[] cf_value = util.process(aux, dt, r, d, q, param_dict, model);
        for (int i = 0; i < cf_value.length; i++) {
            cf_value[i] = cf_value[i].exp();
        }
        Complex[] pF = new Complex[cf_value.length];
        Complex[] dF = new Complex[cf_value.length];
        Complex[] gF = new Complex[cf_value.length];
        for (int i = 0; i < pF.length; i++) {
            pF[i] = cf_value[i].multiply(iutil.multiply(aux[i]).multiply(x - a).exp());
            dF[i] = pF[i].multiply(iutil).multiply(aux[i]);
            gF[i] = pF[i].multiply(iutil.multiply(-1).multiply(aux[i]).add(iutil.multiply(aux[i]).multiply(iutil).multiply(aux[i])));
        }
        pF[0] = pF[0].multiply(0.5);
        dF[0] = dF[0].multiply(0.5);
        gF[0] = gF[0].multiply(0.5);
        double price = 0d, delta = 0d, gamma = 0d;
        for (int i = 0; i < pF.length; i++) {
            price += pF[i].getReal() * V[i][0];
            delta += dF[i].getReal() * V[i][0] / S_0;
            gamma += gF[i].getReal() * V[i][0] / (S_0 * S_0);
        }
        price *= Math.exp(-r * dt);
        delta *= Math.exp(-r * dt);
        gamma *= Math.exp(-r * dt);
        return new double[]{price, delta, gamma};

    }

    protected double[] COS_UO_3(int n, int N_step, double[] c, int cp, double S_0, double HU, double HL, double T, double r, double d,
                                double strike, Map<String, Double> param_dict, int[] UO_list, double mu_f, double P, String model) {
        Complex iutil = new Complex(0, 1);
        double dt = T / N_step;
        double N = Math.pow(2, n);
        double x = Math.log(S_0 / strike);
        double hu = Math.log(HU / strike);
        double hl = Math.log(HL / strike);
        int idx_temp = 0;
//        double cash = (mu_f * T) * P;
        double a = c[0] + x - L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double b = c[0] + x + L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double[][] Grid_i = new double[(int) N][1];
        Complex[] aux = new Complex[(int) N];
        for (int i = 0; i < Grid_i.length; i++) {
            Grid_i[i][0] = i;
            aux[i] = new Complex(pi * i / (b - a));
        }
        if(pro.ID == 12) {
            obs_list = pro.obs_list;
        }
        double[][] t1 = util.Coefficient_V2(Grid_i, new double[]{a}, new double[]{0}, new double[]{a}, new double[]{b});
        double[][] t2 = util.Coefficient_V1(Grid_i, new double[]{0}, new double[]{hu}, new double[]{a}, new double[]{b});
        double[][] V = new double[t1.length][t1[0].length];
        for (int i = 0; i < V.length; i++) {
            for (int j = 0; j < V[0].length; j++) {
                V[i][j] = (t1[i][j] + t2[i][j]) * P;
            }
        }

        if(pro.ID == 1 || pro.ID == 2) {
            double h_limit = Math.log(pro.kd_dict.get("S_T_limit")/strike);
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{0}, new double[]{hu}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V2(Grid_i, new double[]{h_limit}, new double[]{0}, new double[]{a}, new double[]{b});
            double[][] tt3 = util.Coefficient_V1(Grid_i, new double[]{a}, new double[]{h_limit}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = (tt1[i][j] + tt2[i][j] + pro.kd_dict.get("S_T_limit")/S_0 * tt3[i][j]) *P;
                }
            }
        } else if(pro.ID == 3) {
            double hu_2 = Math.log(pro.kd_dict.get("K_2")/strike);
            double hu_1 = Math.log(pro.kd_dict.get("K_1")/strike);
            double[][] tt1 = util.Coefficient_V2(Grid_i, new double[]{hu_1}, new double[]{hu_2}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V1(Grid_i, new double[]{a}, new double[]{hu_1}, new double[]{a}, new double[]{b});
            double[][] tt3 = util.Coefficient_V1(Grid_i, new double[]{hu_1}, new double[]{hu_2}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*P + pro.kd_dict.get("S_T_limit")/S_0*tt2[i][j]*P + (1-pro.kd_dict.get("K_2"))*P*tt3[i][j];
                }
            }
        } else if(pro.ID==5 || pro.ID==6 || pro.ID==7 || pro.ID==13 || pro.ID==15) {
            idx_temp = pro.KO_arr.length - 1;
            HU = pro.KO_arr[idx_temp];
            hu = Math.log(HU/strike);
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{0}, new double[]{hu}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V2(Grid_i, new double[]{a}, new double[]{0}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = (tt1[i][j] + tt2[i][j]) * P;
                }
            }
        } else if(pro.ID == 11) {
            double cash_kd11 = (pro.kd_dict.get("mu_f")*obs_list.length) * P;
            double hk = Math.log(pro.kd_dict.get("K")/strike);
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{hk}, new double[]{b}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V2(Grid_i, new double[]{a}, new double[]{hk}, new double[]{a}, new double[]{b});
            double[][] tt3 = util.Coefficient_V1(Grid_i, new double[]{a}, new double[]{hk}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*(P+cash_kd11) + tt2[i][j]/(pro.kd_dict.get("K")/S_0)*P + tt3[i][j]*cash_kd11;
                }
            }
        } else if(pro.ID == 14) {
            double hk = Math.log(pro.kd_dict.get("K")/strike);
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{hk}, new double[]{hu}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V2(Grid_i, new double[]{a}, new double[]{hk}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j] * P + tt2[i][j]/(pro.kd_dict.get("K")/S_0) * P;
                }
            }
        } else if(pro.ID==16 || pro.ID==18) {
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*(1+pro.kd_dict.get("min_rate")*T - pro.kd_dict.get("opt_val"))*P;
                }
            }
        } else if(pro.ID == 17 || pro.ID == 19) {
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*(1+pro.kd_dict.get("min_rate")*T)*P;
                }
            }
        } else if(pro.ID == 20) {
            double cash_kd20 = (1+mu_f*T)*P;
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V2(Grid_i, new double[]{a}, new double[]{hl}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j] * cash_kd20 + tt2[i][j] * P;
                }
            }
        } else if(pro.ID == 23) {
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{hl}, new double[]{0}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V1(Grid_i, new double[]{0}, new double[]{b}, new double[]{a}, new double[]{b});
            double[][] tt3 = util.Coefficient_V2(Grid_i, new double[]{0}, new double[]{b}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*P + tt2[i][j]*2*P - tt3[i][j]*P;
                }
            }
        }

        for (int m = N_step - 1; m > 0; m--) {
            if (this.contain(obs_list, m)) {
                if (pro.ID==5 || pro.ID==6 || pro.ID==7 || pro.ID==13 || pro.ID==15) {
                    idx_temp -= 1;
                    HU = pro.KO_arr[idx_temp];
                    hu = Math.log(HU/strike);
                    V = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                } else if(pro.ID == 3) {
                    double hu_2 = Math.log(pro.kd_dict.get("K_2")/strike);
                    V = util.C_value(new double[]{a}, new double[]{hu_2}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                } else if(pro.ID == 23) {
                    V = util.C_value(new double[]{hl}, new double[]{b}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                } else {
                    V = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                }
            } else {
                V = util.C_value(new double[]{a}, new double[]{b}, new double[]{a}, new double[]{b},
                        (int) N, V, dt, r, d, q, model, param_dict);
            }
        }
        Complex[] cf_value = util.process(aux, dt, r, d, q, param_dict, model);
        for (int i = 0; i < cf_value.length; i++) {
            cf_value[i] = cf_value[i].exp();
        }
        Complex[] pF = new Complex[cf_value.length];
        Complex[] dF = new Complex[cf_value.length];
        Complex[] gF = new Complex[cf_value.length];
        for (int i = 0; i < pF.length; i++) {
            pF[i] = cf_value[i].multiply(iutil.multiply(aux[i]).multiply(x - a).exp());
            dF[i] = pF[i].multiply(iutil).multiply(aux[i]);
            gF[i] = pF[i].multiply(iutil.multiply(-1).multiply(aux[i]).add(iutil.multiply(aux[i]).multiply(iutil).multiply(aux[i])));
        }
        pF[0] = pF[0].multiply(0.5);
        dF[0] = dF[0].multiply(0.5);
        gF[0] = gF[0].multiply(0.5);
        double price = 0d, delta = 0d, gamma = 0d;
        for (int i = 0; i < pF.length; i++) {
            price += pF[i].getReal() * V[i][0];
            delta += dF[i].getReal() * V[i][0] / S_0;
            gamma += gF[i].getReal() * V[i][0] / (S_0 * S_0);
        }
        price *= Math.exp(-r * dt);
        delta *= Math.exp(-r * dt);
        gamma *= Math.exp(-r * dt);
        if(pro.ID == 10) {
            return new double[]{0,0,0};
        }
        return new double[]{price, delta, gamma};
    }

    protected double[] COS_DUO_4(int n, int N_step, double[] c, int cp, double S_0, double HU, double HL, double T, double r, double d,
                                 double strike, Map<String, Double> param_dict, int[] UO_list, double mu_f, double P, String model) {
        if (pro.ID==2 || pro.ID==3 || pro.ID==11) {
            return new double[]{0,0,0};
        }
        Complex iutil = new Complex(0, 1);
        double dt = T / N_step;
        double N = Math.pow(2, n);
        double x = Math.log(S_0 / strike);
        double hu = Math.log(HU / strike);
        double hl = Math.log(HL / strike);
        //double cash = (mu_f * T) * P;
        int idx_temp = 0;
        double a = c[0] + x - L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double b = c[0] + x + L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double[][] Grid_i = new double[(int) N][1];
        Complex[] aux = new Complex[(int) N];
        for (int i = 0; i < Grid_i.length; i++) {
            Grid_i[i][0] = i;
            aux[i] = new Complex(pi * i / (b - a));
        }
        if (pro.ID == 12) {
            obs_list = pro.obs_list;
        }
        double[][] t1 = util.Coefficient_V2(Grid_i, new double[]{hl}, new double[]{0}, new double[]{a}, new double[]{b});
        double[][] t2 = util.Coefficient_V1(Grid_i, new double[]{0}, new double[]{hu}, new double[]{a}, new double[]{b});
        double[][] V = new double[t1.length][t1[0].length];
        for (int i = 0; i < V.length; i++) {
            for (int j = 0; j < V[0].length; j++) {
                V[i][j] = (t1[i][j] + t2[i][j]) * P;
            }
        }

        if (pro.ID==5 || pro.ID == 6 || pro.ID == 7 || pro.ID==13 || pro.ID == 15) {
            idx_temp = pro.KO_arr.length-1;
            HU = pro.KO_arr[idx_temp];
            hu = Math.log(HU/strike);
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{0}, new double[]{hu}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V2(Grid_i, new double[]{hl}, new double[]{0}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = (tt1[i][j] + tt2[i][j]) * P;
                }
            }
        } else if(pro.ID == 14) {
            double hk = Math.log(pro.kd_dict.get("K")/strike);
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{hk}, new double[]{hu}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V2(Grid_i, new double[]{hl}, new double[]{hk}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*P + tt2[i][j]/(pro.kd_dict.get("K"))*P;
                }
            }
        } else if(pro.ID==16 || pro.ID==18) {
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*(1+pro.kd_dict.get("min_rate")*T - pro.kd_dict.get("opt_val"))*P;
                }
            }
        } else if(pro.ID==17 || pro.ID==19) {
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*(1+pro.kd_dict.get("min_rate")*T)*P;
                }
            }
        } else if(pro.ID==20) {
            double cash_kd20 = (1+mu_f*T)*P;
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*cash_kd20;
                }
            }
        } else if(pro.ID==23) {
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{hl}, new double[]{0}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V1(Grid_i, new double[]{0}, new double[]{hu}, new double[]{a}, new double[]{b});
            double[][] tt3 = util.Coefficient_V2(Grid_i, new double[]{0}, new double[]{hu}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = tt1[i][j]*P + tt2[i][j]*2*P - tt3[i][j]*P;
                }
            }
        } else if(pro.ID==24) {
            HL = pro.KL_arr[pro.KL_arr.length-1];
            hl = Math.log(HL/strike);
            double[][] tt1 = util.Coefficient_V1(Grid_i, new double[]{0}, new double[]{hu}, new double[]{a}, new double[]{b});
            double[][] tt2 = util.Coefficient_V2(Grid_i, new double[]{hl}, new double[]{0}, new double[]{a}, new double[]{b});
            for (int i = 0; i < V.length; i++) {
                for (int j = 0; j < V[0].length; j++) {
                    V[i][j] = (tt1[i][j] + tt2[i][j]) * P;
                }
            }
        }

        for (int m = N_step - 1; m > 0; m--) {
            if (this.contain(obs_list, m)) {
                if(pro.ID == 5 || pro.ID == 6 || pro.ID == 7 || pro.ID == 13 || pro.ID == 15) {
                    idx_temp -= 1;
                    HU = pro.KO_arr[idx_temp];
                    hu = Math.log(HU/strike);
                    V = util.C_value(new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                } else if(pro.ID == 24) {
                    HL = pro.KL_arr[m-1];
                    hl = Math.log(HL/strike);
                    V = util.C_value(new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                } else{
                    V = util.C_value(new double[]{hl}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                }
            } else {
                if(pro.ID == 23) {
                    V = util.C_value(new double[]{a}, new double[]{hu}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                } else if(pro.ID == 24) {
                    HL = pro.KL_arr[m-1];
                    hl = Math.log(HL/strike);
                    V = util.C_value(new double[]{hl}, new double[]{b}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                } else {
                    V = util.C_value(new double[]{hl}, new double[]{b}, new double[]{a}, new double[]{b},
                            (int) N, V, dt, r, d, q, model, param_dict);
                }
            }
        }
        Complex[] cf_value = util.process(aux, dt, r, d, q, param_dict, model);
        for (int i = 0; i < cf_value.length; i++) {
            cf_value[i] = cf_value[i].exp();
        }
        Complex[] pF = new Complex[cf_value.length];
        Complex[] dF = new Complex[cf_value.length];
        Complex[] gF = new Complex[cf_value.length];
        for (int i = 0; i < pF.length; i++) {
            pF[i] = cf_value[i].multiply(iutil.multiply(aux[i]).multiply(x - a).exp());
            dF[i] = pF[i].multiply(iutil).multiply(aux[i]);
            gF[i] = pF[i].multiply(iutil.multiply(-1).multiply(aux[i]).add(iutil.multiply(aux[i]).multiply(iutil).multiply(aux[i])));
        }
        pF[0] = pF[0].multiply(0.5);
        dF[0] = dF[0].multiply(0.5);
        gF[0] = gF[0].multiply(0.5);
        double price = 0d, delta = 0d, gamma = 0d;
        for (int i = 0; i < pF.length; i++) {
            price += pF[i].getReal() * V[i][0];
            delta += dF[i].getReal() * V[i][0] / S_0;
            gamma += gF[i].getReal() * V[i][0] / (S_0 * S_0);
        }
        price *= Math.exp(-r * dt);
        delta *= Math.exp(-r * dt);
        gamma *= Math.exp(-r * dt);

        return new double[]{price, delta, gamma};
    }


    protected double[] COS_kd_25_1(int n, int N_step, double[] c, int cp, double S_0, double HL, double T, double r, double d,
                                   double strike, Map<String, Double> param_dict, int[] UO_list, double mu_f, double P, String model) {
        Complex iutil = new Complex(0, 1);
        double dt = T / N_step;
        double N = Math.pow(2, n);
        double x = Math.log(S_0 / strike);
        double hl = Math.log(HL / strike);
        double a = c[0] + x - L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double b = c[0] + x + L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double[][] Grid_i = new double[(int) N][1];
        Complex[] aux = new Complex[(int) N];
        for (int i = 0; i < Grid_i.length; i++) {
            Grid_i[i][0] = i;
            aux[i] = new Complex(pi * i / (b - a));
        }
        double[][] t1 = util.Coefficient_V1(Grid_i, new double[]{hl}, new double[]{0}, new double[]{a}, new double[]{b});
        double[][] t2 = util.Coefficient_V2(Grid_i, new double[]{0}, new double[]{b}, new double[]{a}, new double[]{b});
        double[][] t3 = util.Coefficient_V1(Grid_i, new double[]{0}, new double[]{b}, new double[]{a}, new double[]{b});
        double[][] V = new double[t1.length][t1[0].length];
        for (int i = 0; i < V.length; i++) {
            for (int j = 0; j < V[0].length; j++) {
                V[i][j] = t1[i][j]*P + (t2[i][j] - t3[i][j])*P*pro.kd_dict.get("NL_rate") + t3[i][j]*P;
            }
        }
        for (int m = N_step - 1; m > 0; m--) {
            V = util.C_value(new double[]{hl}, new double[]{b}, new double[]{a}, new double[]{b},
                    (int) N, V, dt, r, d, q, model, param_dict);
        }
        Complex[] cf_value = util.process(aux, dt, r, d, q, param_dict, model);
        for (int i = 0; i < cf_value.length; i++) {
            cf_value[i] = cf_value[i].exp();
        }
        Complex[] pF = new Complex[cf_value.length];
        Complex[] dF = new Complex[cf_value.length];
        Complex[] gF = new Complex[cf_value.length];
        for (int i = 0; i < pF.length; i++) {
            pF[i] = cf_value[i].multiply(iutil.multiply(aux[i]).multiply(x - a).exp());
            dF[i] = pF[i].multiply(iutil).multiply(aux[i]);
            gF[i] = pF[i].multiply(iutil.multiply(-1).multiply(aux[i]).add(iutil.multiply(aux[i]).multiply(iutil).multiply(aux[i])));
        }
        pF[0] = pF[0].multiply(0.5);
        dF[0] = dF[0].multiply(0.5);
        gF[0] = gF[0].multiply(0.5);
        double price = 0d, delta = 0d, gamma = 0d;
        for (int i = 0; i < pF.length; i++) {
            price += pF[i].getReal() * V[i][0];
            delta += dF[i].getReal() * V[i][0] / S_0;
            gamma += gF[i].getReal() * V[i][0] / (S_0 * S_0);
        }
        price *= Math.exp(-r * dt);
        delta *= Math.exp(-r * dt);
        gamma *= Math.exp(-r * dt);

        return new double[]{price, delta, gamma};

    }

    protected double[] COS_kd_25_2(int n, int N_step, double[] c, int cp, double S_0, double HL, double T, double r, double d,
                                   double strike, Map<String, Double> param_dict, int[] UO_list, double mu_f, double P, String model) {
        Complex iutil = new Complex(0, 1);
        double dt = T / N_step;
        double N = Math.pow(2, n);
        double x = Math.log(S_0 / strike);
        double hl = Math.log(HL / strike);
        double a = c[0] + x - L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double b = c[0] + x + L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double[][] Grid_i = new double[(int) N][1];
        Complex[] aux = new Complex[(int) N];
        for (int i = 0; i < Grid_i.length; i++) {
            Grid_i[i][0] = i;
            aux[i] = new Complex(pi * i / (b - a));
        }
        double[][] t1 = util.Coefficient_V2(Grid_i, new double[]{hl}, new double[]{b}, new double[]{a}, new double[]{b});
        double[][] V = new double[t1.length][t1[0].length];
        for (int i = 0; i < V.length; i++) {
            for (int j = 0; j < V[0].length; j++) {
                V[i][j] = t1[i][j]*P;
            }
        }
        for (int m = N_step - 1; m > 0; m--) {
            V = util.C_value(new double[]{hl}, new double[]{b}, new double[]{a}, new double[]{b},
                    (int) N, V, dt, r, d, q, model, param_dict);
        }
        Complex[] cf_value = util.process(aux, dt, r, d, q, param_dict, model);
        for (int i = 0; i < cf_value.length; i++) {
            cf_value[i] = cf_value[i].exp();
        }
        Complex[] pF = new Complex[cf_value.length];
        Complex[] dF = new Complex[cf_value.length];
        Complex[] gF = new Complex[cf_value.length];
        for (int i = 0; i < pF.length; i++) {
            pF[i] = cf_value[i].multiply(iutil.multiply(aux[i]).multiply(x - a).exp());
            dF[i] = pF[i].multiply(iutil).multiply(aux[i]);
            gF[i] = pF[i].multiply(iutil.multiply(-1).multiply(aux[i]).add(iutil.multiply(aux[i]).multiply(iutil).multiply(aux[i])));
        }
        pF[0] = pF[0].multiply(0.5);
        dF[0] = dF[0].multiply(0.5);
        gF[0] = gF[0].multiply(0.5);
        double price = 0d, delta = 0d, gamma = 0d;
        for (int i = 0; i < pF.length; i++) {
            price += pF[i].getReal() * V[i][0];
            delta += dF[i].getReal() * V[i][0] / S_0;
            gamma += gF[i].getReal() * V[i][0] / (S_0 * S_0);
        }
        price *= Math.exp(-r * dt);
        delta *= Math.exp(-r * dt);
        gamma *= Math.exp(-r * dt);
        return new double[]{price, delta, gamma};
    }

    protected double[] COS_kd_25_3(int n, int N_step, double[] c, int cp, double S_0, double HL, double T, double r, double d,
                                   double strike, Map<String, Double> param_dict, int[] UO_list, double mu_f, double P, String model) {
        Complex iutil = new Complex(0, 1);
        double dt = T / N_step;
        double N = Math.pow(2, n);
        double x = Math.log(S_0 / strike);
        double hl = Math.log(HL / strike);
        double a = c[0] + x - L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double b = c[0] + x + L * Math.sqrt(c[1] + Math.sqrt(c[2]));
        double[][] Grid_i = new double[(int) N][1];
        Complex[] aux = new Complex[(int) N];
        for (int i = 0; i < Grid_i.length; i++) {
            Grid_i[i][0] = i;
            aux[i] = new Complex(pi * i / (b - a));
        }
        double[][] t1 = util.Coefficient_V2(Grid_i, new double[]{a}, new double[]{b}, new double[]{a}, new double[]{b});
        double[][] V = new double[t1.length][t1[0].length];
        for (int i = 0; i < V.length; i++) {
            for (int j = 0; j < V[0].length; j++) {
                V[i][j] = t1[i][j]*P;
            }
        }
        Complex[] cf_value = util.process(aux, T, r, d, q, param_dict, model);
        for (int i = 0; i < cf_value.length; i++) {
            cf_value[i] = cf_value[i].exp();
        }
        Complex[] pF = new Complex[cf_value.length];
        Complex[] dF = new Complex[cf_value.length];
        Complex[] gF = new Complex[cf_value.length];
        for (int i = 0; i < pF.length; i++) {
            pF[i] = cf_value[i].multiply(iutil.multiply(aux[i]).multiply(x - a).exp());
            dF[i] = pF[i].multiply(iutil).multiply(aux[i]);
            gF[i] = pF[i].multiply(iutil.multiply(-1).multiply(aux[i]).add(iutil.multiply(aux[i]).multiply(iutil).multiply(aux[i])));
        }
        pF[0] = pF[0].multiply(0.5);
        dF[0] = dF[0].multiply(0.5);
        gF[0] = gF[0].multiply(0.5);
        double price = 0d, delta = 0d, gamma = 0d;
        for (int i = 0; i < pF.length; i++) {
            price += pF[i].getReal() * V[i][0];
            delta += dF[i].getReal() * V[i][0] / S_0;
            gamma += gF[i].getReal() * V[i][0] / (S_0 * S_0);
        }
        price *= Math.exp(-r * T);
        delta *= Math.exp(-r * T);
        gamma *= Math.exp(-r * T);
        return new double[]{price, delta, gamma};
    }

    public double[] combination_option() {
        if(pro.ID == 25) {
            double[] d25_1 = COS_kd_25_1(n, N_steps, c, 1, S_0, HL, T, r, d, strike, param_dict, null, mu_f, P, model);
            double[] d25_2 = COS_kd_25_2(n, N_steps, c, 1, S_0, HL, T, r, d, strike, param_dict, null, mu_f, P, model);
            double[] d25_3 = COS_kd_25_3(n, N_steps, c, 1, S_0, HL, T, r, d, strike, param_dict, null, mu_f, P, model);
            return new double[]{d25_1[0]+d25_3[0]-d25_2[0], d25_1[1]+d25_3[1]-d25_2[1], d25_1[2]+d25_3[2]-d25_2[2]};
        }
        double[] d1 = this.COS_DUO_1(n, N_steps, c, 1, S_0, HU, HL, T, r, d, strike, param_dict, null, mu_f, P, model);
        double[] d2 = this.COS_UO_2(n, N_steps, c, 1, S_0, HU, T, r, d, strike, param_dict, null, mu_f, P, model);
        double[] d3 = this.COS_UO_3(n, N_steps, c, 1, S_0, HU, HL, T, r, d, strike, param_dict, null, mu_f, P, model);
        double[] d4 = this.COS_DUO_4(n, N_steps, c, 1, S_0, HU, HL, T, r, d, strike, param_dict, null, mu_f, P, model);
        double price = d1[0] + d2[0] + d3[0] - d4[0];
        double delta = d1[1] + d2[1] + d3[1] - d4[1];
        double gamma = d1[2] + d2[2] + d3[2] - d4[2];
        return new double[]{price, delta, gamma};
    }

    protected boolean contain(int[] array, int val) {
        for (int d : array) {
            if (d == val) return true;
        }
        return false;
    }

    private double[][] add(double[][] a1, double[][] a2) {
        double[][] ans = new double[a1.length][a1[0].length];
        for (int i = 0; i < a1.length; i++) {
            for (int j = 0; j < a1[0].length; j++) {
                ans[i][j] = a1[i][j] + a2[i][j];
            }
        }
        return ans;
    }

    public static void main(String[] args) {
        testBS();
        testVG();
        testHeston();
    }

    public static void testBS() {
        FFTCOSStructure classicSnowball = new FFTCOSStructure();
        SNUtil util = new SNUtil();
        classicSnowball.obs_list = new int[]{21, 42, 63, 84, 105, 126, 147, 168, 189, 210, 231, 252};
        classicSnowball.N_steps = 252;
        classicSnowball.T = 1;
        classicSnowball.S_0 = 1;
        classicSnowball.mu_f = 0.2;
        classicSnowball.r = 0.03;
        classicSnowball.q = 0;
        classicSnowball.L = 12;
        classicSnowball.param_dict = new TreeMap<>();
        classicSnowball.param_dict.put("sigma", 0.13);
        classicSnowball.param_dict.put("theta", -0.14);
        classicSnowball.param_dict.put("beta", 0.2);
        classicSnowball.c = util.Cumulants_cal(classicSnowball.r, classicSnowball.q, classicSnowball.T,
                classicSnowball.param_dict, "BS");
        classicSnowball.strike = 1;
        classicSnowball.HU = 1.03;
        classicSnowball.HL = 0.85;
        classicSnowball.n = 10;
        classicSnowball.model = "BS";
        classicSnowball.P = 1;
        double[] res = classicSnowball.combination_option();
        System.out.println("BS: " + Arrays.toString(res));
    }

    public static void testVG() {
        FFTCOSStructure classicSnowball = new FFTCOSStructure();
        SNUtil util = new SNUtil();
        classicSnowball.obs_list = new int[]{21, 42, 63, 84, 105, 126, 147, 168, 189, 210, 231, 252};
        classicSnowball.N_steps = 252;
        classicSnowball.T = 1;
        classicSnowball.S_0 = 1;
        classicSnowball.mu_f = 0.2;
        classicSnowball.r = 0.03;
        classicSnowball.q = 0;
        classicSnowball.L = 12;
        classicSnowball.param_dict = new TreeMap<>();
        classicSnowball.param_dict.put("sigma", 0.13);
        classicSnowball.param_dict.put("theta", -0.14);
        classicSnowball.param_dict.put("beta", 0.2);
        classicSnowball.c = util.Cumulants_cal(classicSnowball.r, classicSnowball.q, classicSnowball.T,
                classicSnowball.param_dict, "VG");
        classicSnowball.strike = 1;
        classicSnowball.HU = 1.03;
        classicSnowball.HL = 0.85;
        classicSnowball.n = 10;
        classicSnowball.model = "VG";
        classicSnowball.P = 1;
        double[] res = classicSnowball.combination_option();
        System.out.println("VG: " + Arrays.toString(res));
    }

    public static void testHeston() {
        FFTCOSStructure classicSnowball = new FFTCOSStructure();
        SNUtil util = new SNUtil();
        classicSnowball.obs_list = new int[]{21, 42, 63, 84, 105, 126, 147, 168, 189, 210, 231, 252};
        classicSnowball.N_steps = 252;
        classicSnowball.T = 1;
        classicSnowball.S_0 = 1;
        classicSnowball.mu_f = 0.2;
        classicSnowball.r = 0.03;
        classicSnowball.q = 0;
        classicSnowball.L = 12;
        classicSnowball.param_dict = new TreeMap<>();
        classicSnowball.param_dict.put("v_0", 0.04);
        classicSnowball.param_dict.put("theta", 0.04);
        classicSnowball.param_dict.put("kappa", 2d);
        classicSnowball.param_dict.put("omega", Math.sqrt(0.02));
        classicSnowball.param_dict.put("rho", -0.7);
        classicSnowball.c = util.Cumulants_cal(classicSnowball.r, classicSnowball.q, classicSnowball.T,
                classicSnowball.param_dict, "Heston");
        classicSnowball.strike = 1;
        classicSnowball.HU = 1.03;
        classicSnowball.HL = 0.85;
        classicSnowball.n = 10;
        classicSnowball.model = "Heston";
        classicSnowball.P = 1;
        double[] res = classicSnowball.combination_option();
        System.out.println("Heston: " + Arrays.toString(res));
    }
}

//0.05085973132144512 -0.29330680481268495 -12.420987920530033
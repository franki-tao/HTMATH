package 结构性存款;

import org.apache.commons.math3.distribution.NormalDistribution;
import 结构性存款.product.product;

import java.util.ArrayList;

/**
 * 基于MC的经典雪球定价
 */
public class MCStructure {
    public double S_0;
    public int N_samples;
    public int N_steps;
    public double T;
    public double r;
    public double q;
    public double sigma;
    public double KI_rate;
    public double KO_rate;
    public double mu_f;
    public int[] obs_list;
    public double P;
    public product pro = new product();
    public final NormalDistribution norm = new NormalDistribution();

    public double[][] GBMpath(int N_samples, int N_steps, double T, double r, double q, double sigma, double S_0) {
        double dt = T / N_steps;
        double[][] S = new double[N_samples][N_steps + 1];
        for (int i = 0; i < N_samples; i++) {
            S[i][0] = S_0;
        }
        double mu = r - q;

        for (int i = 0; i < N_steps; i++) {
            double[] sample = norm.sample(N_samples);
            for (int j = 0; j < N_samples; j++) {
                S[j][i + 1] = S[j][i] + S[j][i] * mu * dt + S[j][i] * sigma * Math.pow(dt, 0.5) * sample[j];
            }
        }
        return S;
    }

    private double V_MC(double S_0, double T, double r, double q, double sigma, int N_samples, int N_steps, double KI_rate, double KO_rate,
                        double mu_f, int[] obs_list, double P) {
        double[][] S_index = GBMpath(N_samples, N_steps, T, r, q, sigma, S_0);
        double[] cf_ko_arr = new double[obs_list.length];
        for (int i = 0; i < cf_ko_arr.length; i++) {
            cf_ko_arr[i] = Math.exp(obs_list[i] / 365d * (-r)) * P * (1 + obs_list[i] / 365d * mu_f);
        }
        if (pro.ID == 12) {
            obs_list = pro.obs_list;
            cf_ko_arr = new double[obs_list.length];
            for (int i = 0; i < cf_ko_arr.length; i++) {
                cf_ko_arr[i] = Math.exp(pro.obs_list[i] / 365d * (-r)) * P * (1 + pro.kd_dict.get("mu_f"));
            }
        } else if (pro.ID == 3) {
            KO_rate = pro.kd_dict.get("K_2");
            for (int i = 0; i < obs_list.length; i++) {
                cf_ko_arr[i] = Math.exp(obs_list[i] / 365d * (-r)) * P * (1 + pro.kd_dict.get("mu_O") * obs_list[i] / 365d);
            }
        } else if (pro.ID == 4 || pro.ID == 7 || pro.ID == 13) {
            for (int i = 0; i < obs_list.length; i++) {
                cf_ko_arr[i] = Math.exp(obs_list[i] / 365d * (-r)) * P * (1 + pro.mu_O[i] * obs_list[i] / 365);
            }
        } else if (pro.ID == 11) {
            for (int i = 0; i < obs_list.length; i++) {
                cf_ko_arr[i] = Math.exp(obs_list[i] / 365d * (-r)) * P * pro.kd_dict.get("mu_f") * (i + 1) + P;
            }
        } else if (pro.ID == 16 || pro.ID == 18) {
            for (int i = 0; i < obs_list.length; i++) {
                cf_ko_arr[i] = Math.exp(obs_list[i] / 365d * (-r)) * P * (1 + obs_list[i] / 365d * pro.kd_dict.get("mu_O") -
                        pro.kd_dict.get("opt_val"));
            }
        } else if (pro.ID == 17 || pro.ID == 19) {
            for (int i = 0; i < obs_list.length; i++) {
                cf_ko_arr[i] = Math.exp(obs_list[i] / 365d * (-r)) * P * (1 + obs_list[i] / 365d * pro.kd_dict.get("mu_O"));
            }
        }

        double[][] S_K0 = new double[N_samples][obs_list.length];

        if (pro.ID == 5 || pro.ID == 6 || pro.ID == 7 || pro.ID == 13 || pro.ID == 15) {
            for (int i = 0; i < N_samples; i++) {
                for (int j = 0; j < obs_list.length; j++) {
                    if (S_index[i][obs_list[j]] >= pro.KO_arr[j]) {
                        S_K0[i][j] = 1;
                    }
                }
            }
        } else if (pro.ID == 23) {
            for (int i = 0; i < N_samples; i++) {
                for (int j = 0; j < obs_list.length; j++) {
                    if (S_index[i][obs_list[j]] <= pro.kd_dict.get("HL")) {
                        S_K0[i][j] = 1;
                    }
                }
            }
        } else {
            for (int i = 0; i < N_samples; i++) {
                for (int j = 0; j < obs_list.length; j++) {
                    if (S_index[i][obs_list[j]] >= KO_rate) {
                        S_K0[i][j] = 1;
                    }
                }
            }
        }


        double[][] S_K1 = new double[N_samples][S_index[0].length];
        if (pro.ID == 23) {
            for (int i = 0; i < N_samples; i++) {
                for (int j = 0; j < S_index[0].length; j++) {
                    if (S_index[i][j] >= KO_rate) {
                        S_K1[i][j] = -2;
                    }
                }
            }
        } else if (pro.ID == 24) {
            for (int i = 0; i < N_samples; i++) {
                for (int j = 0; j < S_index[0].length; j++) {
                    if (S_index[i][j] <= pro.KL_arr[j]) {
                        S_K1[i][j] = -2;
                    }
                }
            }
        } else {
            for (int i = 0; i < N_samples; i++) {
                for (int j = 0; j < S_index[0].length; j++) {
                    if (S_index[i][j] <= KI_rate) {
                        S_K1[i][j] = -2;
                    }
                }
            }
        }


        double[][] S_K0_adjust = new double[N_samples][obs_list.length];
        double[][] S_K1_adjust = new double[N_samples][S_K1[0].length];
        ArrayList<int[]> K0_index = new ArrayList<>();
        //ArrayList<int[]> K1_index = new ArrayList<>();

        for (int i = 0; i < S_K0.length; i++) {
            for (int j = 0; j < S_K0[0].length; j++) {
                if (S_K0[i][j] == 1) {
                    S_K0_adjust[i][j] = 1;
                    K0_index.add(new int[]{i, j});
                    break;
                }
            }
        }

        ArrayList<Double> cf_ko_arr_prt = new ArrayList<>();
        if (pro.ID == 8) {
            for (int[] ko_id : K0_index) {
                cf_ko_arr_prt.add((S_index[ko_id[0]][obs_list[ko_id[1]]] / S_0 - KO_rate / S_0) * pro.kd_dict.get("prt_rate") * P);
            }
        }

        if (pro.ID == 22) {
            double kl_flag = 0d;
            for (int i = 0; i < S_K1.length; i++) {
                kl_flag = 0;
                for (int j = 0; j < S_K1[0].length; j++) {
                    if (S_K1[i][j] == -2) {
                        kl_flag += 1;
                        if (kl_flag == 3) {
                            S_K1_adjust[i][j] = -2;
                            break;
                        }
                    }
                }
            }
        } else {
            for (int i = 0; i < S_K1.length; i++) {
                for (int j = 0; j < S_K1[0].length; j++) {
                    if (S_K1[i][j] == -2) {
                        S_K1_adjust[i][j] = -2;
                        //K1_index.add(new int[]{i, j});
                        break;
                    }
                }
            }
        }


        double[][] S_K0_temp = new double[S_index.length][S_index[0].length];
        for (int i = 0; i < K0_index.size(); i++) {
            int[] temp = K0_index.get(i);
            S_K0_temp[temp[0]][temp[1]] = 1;
        }
        double[][] S_logitic = add(S_K0_temp, S_K1_adjust);

        ArrayList<Double> K1_cf = new ArrayList<>();
        ArrayList<Double> KN_cf = new ArrayList<>();
        ArrayList<Double> NL_cf = new ArrayList<>();
        ArrayList<Double> L_cf = new ArrayList<>();

        for (int i = 0; i < S_logitic.length; i++) {
            double logitic = sum(S_logitic[i]);
            if (pro.ID == 25) {
                if (logitic >= 0) {
                    double cf = Math.exp(-r * T) * (Math.max(S_index[i][S_index[0].length - 1] - S_0, 0) * P *
                            pro.kd_dict.get("NL_rate") + P);
                    NL_cf.add(cf);
                } else {
                    double cf = Math.exp(-r * T) * ((S_index[i][S_index[0].length - 1] - S_0) * P *
                            1 + P);
                    L_cf.add(cf);
                }
            }

            if (logitic == -2) {
                double cf = Math.exp(-r * T) * Math.min(S_index[i][N_steps] / S_0 * P, P);
                if (pro.ID == 1 || pro.ID == 2) {
                    cf = Math.exp(-r * T) * Math.min(Math.max(S_index[i][N_steps], pro.kd_dict.get("S_T_limit")), S_0);
                } else if (pro.ID == 3) {
                    if (pro.kd_dict.get("K_2") > S_index[i][N_steps] && S_index[i][N_steps] > pro.kd_dict.get("K_1")) {
                        cf = Math.exp(-r * T) * (1 - pro.kd_dict.get("K_2") / S_0 + S_index[i][N_steps] / S_0) * P;
                    } else if (pro.kd_dict.get("K_1") > S_index[i][N_steps]) {
                        cf = (pro.kd_dict.get("S_T_limit")) * P;
                    }
                } else if (pro.ID == 11) {
                    cf = Math.min(S_index[i][N_steps] / pro.kd_dict.get("K") * P, P) + P * obs_list.length *
                            pro.kd_dict.get("mu_f");
                } else if (pro.ID == 14) {
                    cf = Math.exp(-r * T) * Math.min(S_index[i][N_steps] / pro.kd_dict.get("K") * P, P);
                } else if (pro.ID == 16 || pro.ID == 18) {
                    cf = Math.exp(-r * T) * P * (1 + pro.kd_dict.get("min_rate") * T - pro.kd_dict.get("opt_val"));
                } else if (pro.ID == 17 || pro.ID == 19) {
                    cf = Math.exp(-r * T) * P * (1 + pro.kd_dict.get("min_rate") * T);
                } else if (pro.ID == 20) {
                    if (S_index[i][N_steps] > KI_rate) {
                        cf = Math.exp(-r * T) * P * (1 + mu_f * T);
                    } else {
                        cf = S_index[i][N_steps] / S_0 * P;
                    }
                } else if (pro.ID == 23) {
                    cf = Math.exp(-r * T) * ((P - Math.max(S_index[i][N_steps] / S_0 * P, P)) + P);
                }
                K1_cf.add(cf);
            } else if (logitic == 0) {
                double cf = Math.exp(-r * T) * P * (1 + mu_f * T);
                if (pro.ID == 2) {
                    cf = Math.exp(-r * T) * Math.min(Math.max(S_index[i][N_steps] / S_0, pro.kd_dict.get("S_T_limit"))
                            / S_0, 1) * P;
                } else if (pro.ID == 3) {
                    if (pro.kd_dict.get("K_2") > S_index[i][N_steps] && S_index[i][N_steps] > pro.kd_dict.get("K_1")) {
                        cf = Math.exp(-r * T) * (1 - pro.kd_dict.get("K_2") / S_0 + S_index[i][N_steps] / S_0) * P;
                    } else if (pro.kd_dict.get("K_1") > S_index[i][N_steps]) {
                        cf = pro.kd_dict.get("S_T_limit") / S_0 * P;
                    }
                } else if (pro.ID == 4 || pro.ID == 7) {
                    cf = Math.exp(-r * T) * P * (1 + pro.kd_dict.get("mu_f"));
                } else if (pro.ID == 9 || pro.ID == 13) {
                    cf = Math.exp(-r * T) * P * (1 + pro.kd_dict.get("mu_f") * T);
                } else if (pro.ID == 11) {
                    cf = Math.min(S_index[i][N_steps] / pro.kd_dict.get("K") * P, P) + P * obs_list.length *
                            pro.kd_dict.get("mu_f");
                } else if (pro.ID == 12) {
                    cf = Math.exp(-r * T) * P * (1 + pro.kd_dict.get("mu_f"));
                } else if (pro.ID == 16) {
                    cf = Math.exp(-r * T) * P * (1 + pro.kd_dict.get("min_rate") * T - pro.kd_dict.get("opt_val"));
                } else if (pro.ID == 17) {
                    cf = Math.exp(-r * T) * P * (1 + pro.kd_dict.get("min_rate") * T);
                } else if (pro.ID == 18) {
                    cf = Math.exp(-r * T) * P * (1 + pro.kd_dict.get("mu_f") * T - pro.kd_dict.get("opt_val"));
                } else if (pro.ID == 19) {
                    cf = Math.exp(-r * T) * P * (1 + pro.kd_dict.get("mu_f") * T);
                }
                KN_cf.add(cf);
            }
        }

        double t3 = 0d;
        for (int i = 0; i < S_K0_adjust.length; i++) {
            double tt = 0d;
            for (int j = 0; j < cf_ko_arr.length; j++) {
                tt += S_K0_adjust[i][j] * cf_ko_arr[j];
            }
            t3 += tt;
        }
        double t1 = 0d;
        for (double d : KN_cf) {
            t1 += d;
        }
        for (double d : K1_cf) {
            t1 += d;
        }
        //t1 += t3;

        if (pro.ID == 8) {
            double t4 = 0;
            for (double d : cf_ko_arr_prt) {
                t4 += d;
            }
            return (t1 + t3 + t4) / N_samples;
        } else if (pro.ID == 25) {
            double t5 = 0d;
            for (double d : L_cf) {
                t5 += d;
            }
            for (double d : NL_cf) {
                t5 += d;
            }
            return (t5) / N_samples;
        }
        return (t1 + t3) / N_samples;
    }

    public double combination_MC() {
        return V_MC(S_0, T, r, q, sigma, N_samples, N_steps, KI_rate, KO_rate, mu_f, obs_list, P);
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

    private double sum(double[] arr) {
        double res = 0d;
        for (double s : arr) {
            res += s;
        }
        return res;
    }

    public static void main(String[] args) {
        MCStructure classicSnowballMC = new MCStructure();
        double v = classicSnowballMC.V_MC(1, 1, 0.03, 0, 0.13, 10000, 252, 0.85, 1.03,
                0.2, new int[]{21, 42, 63, 84, 105, 126, 147, 168, 189, 210, 231, 252}, 1);
        System.out.println(v);
    }
}

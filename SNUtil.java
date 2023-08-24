package 结构性存款;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.special.Gamma;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class SNUtil {
    private final double pi = Math.PI;

    public Complex[] process(Complex[] u, double T, double r, double d, double q, Map<String, Double> param_dict, String model) {
        Complex iunit = new Complex(0, 1); //虚数单位
        Complex[] res = new Complex[u.length];

        if (model.equalsIgnoreCase("BS")) {
            double sigma = param_dict.get("sigma");
            for (int i = 0; i < res.length; i++) {
                res[i] = iunit.multiply(u[i]).multiply(T * (r - d - 0.5 * sigma * sigma)).
                        subtract(u[i].multiply(u[i]).multiply(0.5 * sigma * sigma * T));
            }
        } else if (model.equalsIgnoreCase("Heston")) {
            double v_0 = param_dict.get("v_0");
            double theta = param_dict.get("theta");
            double kappa = param_dict.get("kappa");
            double omega = param_dict.get("omega");
            double rho = param_dict.get("rho");
            double gamma = 0.5 * omega * omega;
            for (int i = 0; i < res.length; i++) {
                Complex alpha = (u[i].multiply(u[i]).add(u[i].multiply(iunit))).multiply(-0.5);
                Complex beta = u[i].multiply(iunit).multiply(-rho * omega).add(kappa);
                Complex D = beta.multiply(beta).subtract(alpha.multiply(4 * gamma)).sqrt();
                Complex bD = beta.subtract(D);
                Complex expDT = D.multiply(-T).exp();
                Complex G = bD.divide(beta.add(D));
                Complex B = (beta.subtract(D)).divide(omega * omega).multiply(expDT.multiply(-1).add(1).divide(G.multiply(expDT.multiply(-1)).add(1)));
                Complex psi = G.multiply(expDT).subtract(1).divide(G.subtract(1));
                Complex A = new Complex(kappa * theta / omega / omega).multiply(beta.subtract(D).multiply(T).subtract(psi.log().multiply(2)));
                res[i] = A.add(B.multiply(v_0)).add(iunit.multiply(u[i]).multiply((r - d) * T));
            }
        } else if (model.equalsIgnoreCase("VG")) {
            double sigma = param_dict.get("sigma");
            double beta = param_dict.get("beta");
            double theta = param_dict.get("theta");
            double omega = (1 / beta) * Math.log(1 - theta * beta - sigma * sigma * beta * 0.5);
            for (int i = 0; i < res.length; i++) {
                Complex temp = u[i].multiply(iunit).multiply(-theta * beta).add(1).add(u[i].multiply(u[i]).multiply(0.5 * sigma * sigma * beta));
                res[i] = iunit.multiply(u[i]).multiply(T * (r + omega - d)).subtract(new Complex(T).divide(beta).multiply(temp.log()));
            }
        } else if (model.equalsIgnoreCase("CGMY")) {
            double c = param_dict.get("C");
            double g = param_dict.get("G");
            double m = param_dict.get("M");
            double y = param_dict.get("Y");
            double sigma = param_dict.get("sigma");
            double omega = -c * Gamma.gamma(-y) * (Math.pow(m - 1, y) - Math.pow(m, y) + Math.pow(g + 1, y) - Math.pow(g, y));
            for (int i = 0; i < res.length; i++) {
                Complex temp = new Complex(c * T * Gamma.gamma(-y)).multiply(new Complex(m).subtract(iunit.multiply(u[i])).pow(y).subtract(Math.pow(m, y)).add(
                        new Complex(g).add(iunit.multiply(u[i])).pow(y)).subtract(Math.pow(g, y)));
                res[i] = iunit.multiply(u[i]).multiply(T * (r - d + omega - 0.5 * sigma * sigma)).subtract(u[i].multiply(u[i]).
                        multiply(0.5 * sigma * sigma * T)).add(temp);
            }
        } else if (model.equalsIgnoreCase("NIG")) {
            double alpha = param_dict.get("alpha");
            double beta = param_dict.get("beta");
            double delta = param_dict.get("delta");
            double sigma = param_dict.get("sigma");
            double omega = delta * (Math.sqrt(alpha * alpha - (beta + 1) * (beta + 1)) - Math.sqrt(alpha * alpha - beta * beta));
            double mu = r - 0.5 * sigma * sigma - q + omega;
            for (int i = 0; i < res.length; i++) {
                Complex temp = new Complex(alpha * alpha - beta * beta).sqrt().subtract(new Complex(alpha * alpha).
                        subtract(iunit.multiply(u[i]).add(beta).pow(2)).sqrt()).multiply(delta * T);
                res[i] = iunit.multiply(u[i]).multiply(mu * T).subtract(u[i].multiply(u[i]).multiply(0.5 * sigma * sigma * T)).add(temp);
            }
        }
        return res;
    }

    public double[] Cumulants_cal(double r, double q, double tau, Map<String, Double> param_dict, String model) {
        double c1 = 0, c2 = 0, c3 = 0;
        if (model.equalsIgnoreCase("BS")) {
            double sigma = param_dict.get("sigma");
            c1 = (r - q - 0.5 * sigma * sigma) * tau;
            c2 = sigma * sigma * tau;
        } else if (model.equalsIgnoreCase("Heston")) {
            double v_0 = param_dict.get("v_0");
            double theta = param_dict.get("theta");
            double kappa = param_dict.get("kappa");
            double omega = param_dict.get("omega");
            double rho = param_dict.get("rho");
            double gamma = 0.5 * omega * omega;
            c1 = (r - q) * tau + (1 - Math.exp(-kappa * tau)) * (theta - v_0) / (2 * kappa) - 0.5 * theta * tau;
            double term1 = gamma * tau * kappa * Math.exp(-kappa * tau) * (v_0 - theta) * (8 * kappa * rho - 4 * gamma);
            double term2 = kappa * rho * gamma * (1 - Math.exp(-kappa * tau)) * (16 * theta - 8 * v_0);
            double term3 = 2 * theta * kappa * tau * (-4 * kappa * rho * gamma + gamma * gamma + 4 * kappa * kappa);
            double term4 = gamma * gamma * ((theta - 2 * v_0) * Math.exp(-2 * kappa * tau) + theta * (6 * Math.exp(-kappa * tau) - 7) + 2 * v_0);
            double term5 = 8 * kappa * kappa * (v_0 - theta) * (1 - Math.exp(-kappa * tau));
            c2 = 1 / (8 * kappa * kappa * kappa) * (term1 + term2 + term3 + term4 + term5);

        } else if (model.equalsIgnoreCase("VG")) {
            double sigma = param_dict.get("sigma");
            double beta = param_dict.get("beta");
            double theta = param_dict.get("theta");
            double omega = (1 / beta) * Math.log(1 - theta * beta - sigma * sigma * beta * 0.5);

            c1 = (r - q - omega + theta) * tau;
            c2 = (sigma * sigma + beta * theta * theta) * tau;
            c3 = 3 * (sigma * sigma * beta + 2 * Math.pow(theta, 4) * beta * beta * beta +
                    4 * sigma * sigma * theta * theta * beta * beta) * tau;
        } else if (model.equalsIgnoreCase("CGMY")) {
            double C = param_dict.get("C");
            double G = param_dict.get("G");
            double M = param_dict.get("M");
            double Y = param_dict.get("Y");
            double sigma = param_dict.get("sigma");
            double omega = -C * Gamma.gamma(-Y) * (Math.pow(M - 1, Y) - Math.pow(M, Y) + Math.pow(G + 1, Y) - Math.pow(G, Y));

            c1 = (r - q + omega - 0.5 * sigma * sigma) * tau + C * tau * Gamma.gamma(1 - Y) * (Math.pow(M, (Y - 1)) - Math.pow(G, (Y - 1)));
            c2 = sigma * sigma * tau + C * tau * Gamma.gamma(2 - Y) * (Math.pow(M, (Y - 2)) + Math.pow(G, (Y - 2)));
            c3 = C * Gamma.gamma(4 - Y) * tau * (Math.pow(M, (Y - 4)) + Math.pow(G, (Y - 4)));
        } else if (model.equalsIgnoreCase("NIG")) {
            double alpha = param_dict.get("alpha");
            double beta = param_dict.get("beta");
            double delta = param_dict.get("delta");
            double sigma = param_dict.get("sigma");
            double omega = delta * (Math.sqrt(alpha * alpha - (beta + 1) * (beta + 1)) - Math.sqrt(alpha * alpha - beta * beta));
            c1 = (r - q + omega - 0.5 * sigma * sigma + delta * beta / Math.sqrt(alpha * alpha - beta * beta)) * tau;
            c2 = delta * alpha * alpha * tau * Math.pow(alpha * alpha - beta * beta, -3d / 2);
            c3 = 3 * delta * alpha * alpha * (alpha * alpha + 4 * beta * beta) * tau * Math.pow(alpha * alpha - beta * beta, -7d / 2);
        }
        return new double[]{c1, c2, c3};
    }

    public ArrayList<double[][]> Coefficient_b(double[][] k, double[] x_1, double[] x_2, double[] a, double[] b) {
        int N = k.length;
        int M = x_1.length;
        double[][] arg2 = new double[N][M];
        double[][] arg1 = new double[N][M];
        double[][] temp1 = new double[N][M];
        double[][] temp2 = new double[N][M];
        double[][] temp3 = new double[N][M];
        double[][] temp4 = new double[N][M];
        double[][] chi = new double[N][M];
        double[][] psi = new double[N][M];

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                arg2[i][j] = k[i][j] * pi * (x_2[j] - a[j]) / (b[j] - a[j]);
                arg1[i][j] = k[i][j] * pi * (x_1[j] - a[j]) / (b[j] - a[j]);
                temp1[i][j] = Math.cos(arg2[i][j]) * Math.exp(x_2[j]);
                temp2[i][j] = Math.cos(arg1[i][j]) * Math.exp(x_1[j]);
                temp3[i][j] = pi * k[i][j] * Math.sin(arg2[i][j]) * Math.exp(x_2[j]) / (b[j] - a[j]);
                temp4[i][j] = pi * k[i][j] * Math.sin(arg1[i][j]) * Math.exp(x_1[j]) / (b[j] - a[j]);
                chi[i][j] = 1d / (1 + Math.pow((k[i][j] * pi / (b[j] - a[j])), 2)) *
                        (temp1[i][j] - temp2[i][j] + temp3[i][j] - temp4[i][j]);
                psi[i][j] = (Math.sin(arg2[i][j]) - Math.sin(arg1[i][j])) / (k[i][j] * pi) * (b[j] - a[j]);
            }

        }
        for (int i = 0; i < M; i++) {
            chi[0][i] = Math.exp(x_2[i]) - Math.exp(x_1[i]);
            psi[0][i] = x_2[i] - x_1[i];
        }

        ArrayList<double[][]> ans = new ArrayList<>();
        ans.add(chi);
        ans.add(psi);
        return ans;
    }

    public double[][] C_value(double[] x_1, double[] x_2, double[] a, double[] b, int N, double[][] V, double t, double r, double d,
                              double q, String model, Map<String, Double> param_dict) {
        int Nstrike = x_1.length;
        double sigma;
        if (model.equalsIgnoreCase("BS")) {
            sigma = param_dict.get("sigma");
        }
        double[][] temp = new double[N][Nstrike];
        for (int i = 1; i < N + 1; i++) {
            for (int j = 0; j < Nstrike; j++) {
                temp[i - 1][j] = i;
            }
        }
        Complex ii = new Complex(0, 1);
        Complex[][] exp2 = new Complex[N][Nstrike];
        Complex[][] exp1 = new Complex[N][Nstrike];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < Nstrike; j++) {
                exp2[i][j] = ii.multiply(temp[i][j] * pi * (x_2[j] - a[j]) / (b[j] - a[j])).exp();
                exp1[i][j] = ii.multiply(temp[i][j] * pi * (x_1[j] - a[j]) / (b[j] - a[j])).exp();
            }
        }
        Complex[][] m = new Complex[3 * N - 1][Nstrike];
        //N-1
        for (int i = 0; i < Nstrike; i++) {
            m[N - 1][i] = ii.multiply(pi * (x_2[i] - x_1[i]) / (b[i] - a[i]));
        }
        //N~2N
        for (int i = N; i < 2 * N; i++) {
            for (int j = 0; j < Nstrike; j++) {
                m[i][j] = exp2[i - N][j].subtract(exp1[i - N][j]).divide(temp[i - N][j]);
                //exp2[i - N][j].subtract(exp1[i - N][j]).multiply(temp[i - N][j])
            }
        }
        //0~N-1
        for (int i = 0; i < N - 1; i++) {
            for (int j = 0; j < Nstrike; j++) {
                m[i][j] = m[2 * N - 2 - i][j].conjugate().negate();
            }
        }
        //2*N~3*N-1
        for (int i = 2 * N; i < 3 * N - 1; i++) {
            for (int j = 0; j < Nstrike; j++) {
                m[i][j] = exp2[i - 2 * N][j].multiply(exp2[N - 1][j]).subtract(exp1[i - 2 * N][j].
                        multiply(exp1[N - 1][j])).divide(temp[i - 2 * N][j] + N);
            }
        }
        double[][] u_cf = new double[N][Nstrike];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < Nstrike; j++) {
                u_cf[i][j] = (temp[i][j] - 1) * pi / (b[j] - a[j]);
            }
        }

        Complex[] u_cf0 = new Complex[N];
        for (int i = 0; i < N; i++) {
            u_cf0[i] = new Complex(u_cf[i][0]);
        }
        Complex[] process = this.process(u_cf0, t, r, d, q, param_dict, model);

        Complex[][] u = new Complex[N][Nstrike];

        for (int i = 0; i < N; i++) {
            u[i][0] = process[i].exp().multiply(V[i][0]);
        }

        for (int i = 0; i < Nstrike; i++) {
            u[0][i] = u[0][i].multiply(0.5d);
        }


        Complex[][] m_s = new Complex[2 * N][Nstrike];
        //0~N
        for (int i = N - 1; i > -1; i--) {
            System.arraycopy(m[i], 0, m_s[-i + N - 1], 0, Nstrike);
        }
        //N
        for (int i = 0; i < Nstrike; i++) {
            m_s[N][i] = new Complex(0, 0);
        }
        //N+1~2N+1
        for (int i = 2 * N - 2; i > N - 1; i--) {
            System.arraycopy(m[i], 0, m_s[-i + 3 * N - 1], 0, Nstrike);
        }

        Complex[][] u_s = new Complex[2 * N][Nstrike];
        for (int i = 0; i < 2 * N; i++) {
            for (int j = 0; j < Nstrike; j++) {
                if (i < N) {
                    u_s[i][j] = u[i][j];
                } else {
                    u_s[i][j] = new Complex(0, 0);
                }
            }
        }

        Complex[][] m_c = new Complex[2 * N][Nstrike];
        for (int i = 0; i < 2 * N; i++) {
            System.arraycopy(m[3 * N - 2 - i], 0, m_c[i], 0, Nstrike);
        }

        double[][] zeta = new double[2 * N][Nstrike];
        for (int i = 0; i < 2 * N; i++) {
            for (int j = 0; j < Nstrike; j++) {
                zeta[i][j] = i % 2 == 0 ? 1 : -1;
            }
        }

        Complex[] fft_xi_s = new Complex[2 * N];
        Complex[] fft_xi_c = new Complex[2 * N];


        FFT fft = new FFT(getline(u_s, 0));
        Complex[] fft_u_s = fft.fft();
        FFT fft1 = new FFT(getline(m_s, 0));
        Complex[] ht1 = fft1.fft();
        for (int i = 0; i < 2 * N; i++) {
            fft_xi_s[i] = ht1[i].multiply(fft_u_s[i]);
        }
        FFT fft2 = new FFT(fft_xi_s);
        Complex[] xi_s = fft2.ifft();

        FFT fft3 = new FFT(getline(m_c, 0));
        Complex[] ht3 = fft3.fft();

        for (int i = 0; i < 2 * N; i++) {
            fft_xi_c[i] = ht3[i].multiply(fft_u_s[i].multiply(zeta[i][0]));
        }
        FFT fft4 = new FFT(fft_xi_c);
        Complex[] xi_c = fft4.ifft();

        double[][] ans = new double[N][Nstrike];

        for (int i = 0; i < N; i++) {
//            System.out.println(xi_c[i]);
            for (int j = 0; j < Nstrike; j++) {
                ans[i][j] = Math.exp(-r * t) / pi * ((xi_s[i].add(xi_c[N - i - 1])).getImaginary());
            }
        }
        return ans;
    }

    public double[][] Coefficient_V(double[][] k, double[] x_1, double[] x_2, double[] a, double[] b, double cash) {
        ArrayList<double[][]> re = this.Coefficient_b(k, x_1, x_2, a, b);
        double[][] psi = re.get(1);
        double[][] C = new double[a.length][a.length];
        for (int i = 0; i < C.length; i++) {
            C[i][i] = cash / (b[i] - a[i]);
        }
        double[][] result = new double[psi.length][C.length];
        for (int i = 0; i < psi.length; i++) {
            for (int j = 0; j < C.length; j++) {
                result[i][j] = 2 * psi[i][j] * C[j][j];
            }
        }
        return result;
    }

    public double[][] Coefficient_V1(double[][] k, double[] x_1, double[] x_2, double[] a, double[] b) {
        ArrayList<double[][]> re = this.Coefficient_b(k, x_1, x_2, a, b);
        double[][] psi = re.get(1);
        double[][] C = new double[a.length][a.length];
        for (int i = 0; i < C.length; i++) {
            C[i][i] = 1d / (b[i] - a[i]);
        }
        double[][] result = new double[psi.length][C.length];
        for (int i = 0; i < psi.length; i++) {
            for (int j = 0; j < C.length; j++) {
                result[i][j] = 2 * psi[i][j] * C[j][j];
            }
        }
        return result;
    }

    public double[][] Coefficient_V2(double[][] k, double[] x_1, double[] x_2, double[] a, double[] b) {
        ArrayList<double[][]> re = this.Coefficient_b(k, x_1, x_2, a, b);
        double[][] chi = re.get(0);
        double[][] C = new double[a.length][a.length];
        for (int i = 0; i < C.length; i++) {
            C[i][i] = 1d / (b[i] - a[i]);
        }
        double[][] result = new double[chi.length][C.length];
        for (int i = 0; i < chi.length; i++) {
            for (int j = 0; j < C.length; j++) {
                result[i][j] = 2 * chi[i][j] * C[j][j];
            }
        }
        return result;
    }


    private Complex[] getline(Complex[][] data, int index) {
        Complex[] ans = new Complex[data.length];
        for (int i = 0; i < data.length; i++) {
            ans[i] = data[i][index];
        }
        return ans;
    }

    public double[][] Coefficient_V_kd8(double[][] k, double[] x_1, double[] x_2, double[] a, double[] b) {
        ArrayList<double[][]> doubles = Coefficient_b(k, x_1, x_2, a, b);
        double[][] chi = doubles.get(0);
        double[][] ans = new double[chi.length][a.length];
        for (int i = 0; i < chi.length; i++) {
            for (int j = 0; j < a.length; j++) {
                ans[i][j] = 2 * chi[i][j] / (b[j] - a[j]);
            }
        }
        return ans;
    }
}

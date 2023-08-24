package 结构性存款.product;

import org.apache.hadoop.yarn.webapp.hamlet2.Hamlet;
import 结构性存款.FFTCOSStructure;
import 结构性存款.SNUtil;

import java.util.Arrays;
import java.util.TreeMap;

//Out-of-the-Money Call Limited Loss Snowball
public class OMCLLSnowBall {
    public final String product_Name = "虚值看涨限亏雪球";
    public final int ID = 3;
    public product pro;
    public double S_T_limit;
    public double K_2; //行权价2
    public double K_1; //行权价1
    public double mu_O; //敲出票息

    public OMCLLSnowBall() {
    }

    public OMCLLSnowBall(double s_T_limit, double k_2, double k_1, double mu_O) {
        S_T_limit = s_T_limit;
        K_2 = k_2;
        K_1 = k_1;
        this.mu_O = mu_O;
        pro = new product();
        pro.ID = this.ID;
        pro.kd_dict = new TreeMap<>();
        pro.kd_dict.put("S_T_limit", s_T_limit);
        pro.kd_dict.put("K_2", k_2);
        pro.kd_dict.put("K_1", k_1);
        pro.kd_dict.put("mu_O", mu_O);
    }

    public double[] cal(FFTCOSStructure structure) {
        structure.pro = this.pro;
        return structure.combination_option();
    }

    //测试
    public static void main(String[] args) {
        FFTCOSStructure pro = new FFTCOSStructure();
        SNUtil util = new SNUtil();
        pro.obs_list = new int[]{31, 61, 92, 122, 152, 183, 213, 244, 274, 304, 335, 365};
        pro.N_steps = 365;
        pro.T = 1;
        pro.S_0 = 1;
        pro.mu_f = 0.21;
        pro.r = 0.021916;
        pro.q = 0.054633;
        pro.L = 12;
        pro.param_dict = new TreeMap<>();
        pro.param_dict.put("sigma", 0.2202);
        pro.c = util.Cumulants_cal(pro.r, pro.q, pro.T,
                pro.param_dict, "BS");
        pro.strike = 1;
        pro.HU = 1.03;
        pro.HL = 0.8;
        pro.n = 10;
        pro.model = "BS";
        pro.P = 1;
        pro.d = 0.054633;

        OMCLLSnowBall snowBall = new OMCLLSnowBall(0.97, 1.05, 1.02, 0.0195);
        double[] cal = snowBall.cal(pro);
        System.out.println(Arrays.toString(cal)); //[0.9890074402991371, 0.2537407345033089, 0.17710941723562845]
    }
}

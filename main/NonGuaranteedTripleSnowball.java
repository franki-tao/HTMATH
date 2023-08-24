package 结构性存款.product;

import 结构性存款.FFTCOSStructure;
import 结构性存款.SNUtil;

import java.util.Arrays;
import java.util.TreeMap;

public class NonGuaranteedTripleSnowball {
    public final String product_Name = "非保本三元雪球";
    public final int ID = 18;
    public product pro;
    public double min_rate;
    public double mu_O;
    public double mu_f;
    public double opt_val;

    public NonGuaranteedTripleSnowball(double min_rate, double mu_O, double mu_f, double opt_val) {
        this.min_rate = min_rate;
        this.mu_O = mu_O;
        this.mu_f = mu_f;
        this.opt_val = opt_val;
        pro = new product();
        pro.ID = ID;
        pro.kd_dict = new TreeMap<>();
        pro.kd_dict.put("min_rate", min_rate);
        pro.kd_dict.put("mu_O", mu_O);
        pro.kd_dict.put("mu_f", mu_f);
        pro.kd_dict.put("opt_val", opt_val);
    }

    public double[] cal(FFTCOSStructure structure) {
        structure.pro = pro;
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

        NonGuaranteedTripleSnowball snowball = new NonGuaranteedTripleSnowball(
                0.001,
                0.06,
                0.04,
                0.01
        );
        System.out.println(Arrays.toString(snowball.cal(pro))); //[0.9924261990709318, 0.031628733348498916, -0.6489396171818673]
    }

}

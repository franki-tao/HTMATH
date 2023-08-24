package 结构性存款.product;

import 结构性存款.FFTCOSStructure;
import 结构性存款.SNUtil;

import java.util.Arrays;
import java.util.TreeMap;

public class TriggeredSnowball {
    public final String product_Name = "触发型雪球";
    public final int ID = 12;
    public product pro;
    public int[] obs_list;
    public double mu_f;

    public TriggeredSnowball(int[] obs_list, double mu_f) {
        this.obs_list = obs_list;
        this.mu_f = mu_f;
        pro = new product();
        pro.ID = this.ID;
        pro.obs_list = obs_list;
        pro.kd_dict = new TreeMap<>();
        pro.kd_dict.put("mu_f", mu_f);
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

        TriggeredSnowball snowball = new TriggeredSnowball(
                new int[]{7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 112, 119,
                        126, 133, 140, 147, 154, 161, 168, 175, 182, 189, 196, 203, 210, 217, 224,
                        231, 238, 245, 252, 259, 266, 273, 280, 287, 294, 301, 308, 315, 322, 329, 336,
                        343, 350, 357, 364},
                0.08
        );
        System.out.println(Arrays.toString(snowball.cal(pro))); //[1.0130498880656167, 1.2812041460754, -3.0385932847322352]
    }

}

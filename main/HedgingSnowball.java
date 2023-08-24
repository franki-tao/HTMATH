package 结构性存款.product;

import 结构性存款.FFTCOSStructure;
import 结构性存款.SNUtil;

import java.util.Arrays;
import java.util.TreeMap;

public class HedgingSnowball {
    public final String product_Name = "避险雪球";
    public final int ID = 15;
    public product pro;
    public double[] KO_arr;

    public HedgingSnowball(double[] KO_arr) {
        this.KO_arr = KO_arr;
        pro = new product();
        pro.ID = ID;
        pro.KO_arr = KO_arr;
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

        HedgingSnowball snowball = new HedgingSnowball(
                new double[]{1.03, 1.03, 1.03, 1.03, 1.03, 0.8, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03}
        );
        System.out.println(Arrays.toString(snowball.cal(pro))); //[1.0217614369737174, 0.11544612551411049, -5.2807422138003774]
    }

}

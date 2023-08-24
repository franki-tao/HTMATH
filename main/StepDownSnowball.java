package 结构性存款.product;

import 结构性存款.FFTCOSStructure;
import 结构性存款.SNUtil;

import java.util.Arrays;
import java.util.TreeMap;

public class StepDownSnowball {
    public final String product_Name = "降敲雪球";
    public final int ID = 5;
    public product pro;
    public double[] KO_arr; //敲出价列表

    public StepDownSnowball(double[] KO_arr) {
        this.KO_arr = KO_arr;
        pro = new product();
        pro.ID = this.ID;
        pro.KO_arr = KO_arr;
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

        StepDownSnowball snowball = new StepDownSnowball(new double[]{1.03, 1.02, 1.01, 1, 0.99, 0.98, 0.97, 0.96, 0.95,
                0.94, 0.93, 0.92});
        System.out.println(Arrays.toString(snowball.cal(pro))); //[0.9876089498009869, 0.47834922938007485, -5.488935210633542]
    }
}

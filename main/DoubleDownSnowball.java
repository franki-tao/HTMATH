package 结构性存款.product;

import 结构性存款.FFTCOSStructure;
import 结构性存款.SNUtil;

import java.util.Arrays;
import java.util.TreeMap;

public class DoubleDownSnowball {
    public final String product_Name = "双降雪球";
    public final int ID = 7;
    public product pro;
    public double[] KO_arr; //敲出价列表
    public double[] mu_O;
    public double mu_f;

    public DoubleDownSnowball(double[] KO_arr, double[] mu_O, double mu_f) {
        this.KO_arr = KO_arr;
        this.mu_O = mu_O;
        this.mu_f = mu_f;
        pro = new product();
        pro.ID = this.ID;
        pro.kd_dict = new TreeMap<>();
        pro.kd_dict.put("mu_f", mu_f);
        pro.KO_arr = KO_arr;
        pro.mu_O = mu_O;
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
        DoubleDownSnowball snowball = new DoubleDownSnowball(
                new double[]{1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 0.95},
                new double[]{0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08},
                0.08
        );
        System.out.println(Arrays.toString(snowball.cal(pro))); //[0.9684855735953413, 0.7000750821155913, -5.363213012682883]
    }
}

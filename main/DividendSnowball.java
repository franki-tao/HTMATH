package 结构性存款.product;

import 结构性存款.FFTCOSStructure;
import 结构性存款.SNUtil;

import java.util.Arrays;
import java.util.TreeMap;

public class DividendSnowball {
    public final String product_Name = "红利雪球";
    public final int ID = 9;
    public product pro;
    public double mu_f;

    public DividendSnowball(double mu_f) {
        this.mu_f = mu_f;
        pro = new product();
        pro.ID = this.ID;
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

        DividendSnowball snowball = new DividendSnowball(0.08);
        System.out.println(Arrays.toString(snowball.cal(pro))); //[0.9724313243925968, 0.6334307413900641, -5.4379588880452125]
    }
}

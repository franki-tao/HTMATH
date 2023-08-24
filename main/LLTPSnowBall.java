package 结构性存款.product;

import 结构性存款.FFTCOSStructure;
import 结构性存款.SNUtil;

import java.util.Arrays;
import java.util.TreeMap;

public class LLTPSnowBall {
    public final String product_Name = "限亏止盈雪球";
    public final int ID = 2;
    public product pro;
    public double S_T_limit; //保底水平

    public LLTPSnowBall() {
    }

    public LLTPSnowBall(double S_T_limit) {
        this.S_T_limit = S_T_limit;
        pro = new product();
        pro.ID = this.ID;
        pro.kd_dict = new TreeMap<>();
        pro.kd_dict.put("S_T_limit", S_T_limit);
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

        LLTPSnowBall snowball = new LLTPSnowBall(0.8);
        double[] cal = snowball.cal(pro);
        System.out.println(Arrays.toString(cal)); //[0.9781389048440046, 0.5196082711295742, -3.6523465334441507]
    }
}

package 结构性存款.product;

import org.apache.hadoop.yarn.webapp.hamlet2.Hamlet;
import 结构性存款.FFTCOSStructure;
import 结构性存款.SNUtil;

import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

//经典雪球
public class ClassicSnowBall {
    public final String product_Name = "经典雪球";
    public final int ID = 0;
    public product pro;
    public ClassicSnowBall() {
        pro = new product();
        pro.ID = ID;
    }
    public double[] cal(FFTCOSStructure structure) {
        structure.pro = this.pro;
        return structure.combination_option();
    }

    //測試經典雪球
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

        ClassicSnowBall snowBall = new ClassicSnowBall();
        double[] cal = snowBall.cal(pro);
        System.out.println(Arrays.toString(cal)); //[0.9799967064842405, 0.568139290692299, -5.855124940168745]
    }
}

package ht_Math.Complex;
/*
作为数学的基本运算工具 复数
此类包含复数的常见操作， 包括 +， - * /
常见的函数操作， 取指数， 对数， 三角函数，共轭， 取模等
计算基于欧拉公式 a+bi = r*exp(i*theta) = r*(cos(theta) + i*sin(theta))
其中 i为虚数单位， r为模长， theta为辐角 取arctan(b, a)；
 */


public class complex {
    private double Re; //实部
    private double Im; //虚部
    private double FX = 1e-10; //容差

    public complex() {
        Re = 0;
        Im = 0;
    }

    public complex(double re, double im) {
        Re = re;
        Im = im;
    }

    public complex add(complex cmp) {
        return new complex(this.Re + cmp.Re, this.Im + cmp.Im);
    }

    public complex add(double hr) {
        return new complex(this.Re + hr, this.Im);
    }

    public complex subtract(complex cmp) {
        return new complex(this.Re - cmp.Re, this.Im - cmp.Im);
    }

    public complex subtract(double hr) {
        return new complex(this.Re - hr, this.Im);
    }

    public complex multiply(complex cmp) {
        return new complex(this.Re * cmp.Re - this.Im * cmp.Im, this.Re * cmp.Im + cmp.Re * this.Im);
    }

    public complex multiply(double hr) {
        return new complex(this.Re * hr, this.Im * hr);
    }

    public complex division(complex cmp) {
        double tmp = cmp.Re * cmp.Re + cmp.Im * cmp.Im;
        return new complex((this.Re * cmp.Re + this.Im * cmp.Im) / tmp, (-this.Re * cmp.Im + cmp.Re * this.Im) / tmp);
    }

    public complex division(double hr) {
        return new complex(this.Re / hr, this.Im / hr);
    }

    public double Mod() {
        return Math.sqrt(Re * Re + Im * Im);
    }

    public complex conjugate() {
        return new complex(this.Re, -this.Im);
    }

    public complex exp() {
        return new complex(Math.exp(Re) * Math.cos(Im), Math.exp(Re) * Math.sin(Im));
    }

    public complex pow(double index) {
        double tp1 = Mod();
        double arg = Math.atan2(Im, Re);
        double hre = Math.pow(tp1, index) * Math.cos(index * arg);
        double him = Math.pow(tp1, index) * Math.sin(index * arg);
        double floorRe = Math.floor(hre), floorIm = Math.floor(him);
        double ceilRe = Math.ceil(hre), ceilIm = Math.ceil(him);
        if (Math.abs(hre - floorRe) < FX) {
            hre = floorRe;
        }

        if (Math.abs(hre - ceilRe) < FX) {
            hre = ceilRe;
        }

        if (Math.abs(him - floorIm) < FX) {
            him = floorIm;
        }

        if (Math.abs(him - ceilIm) < FX) {
            him = ceilIm;
        }

        return new complex(hre, him);
    }

    public complex sqrt() {
        return pow(0.5);
    }

    public complex negative() {
        return new complex(-Re, -Im);
    }

    public complex sin() {
//        实部 = sin(a) * cosh(b)
//        虚部 = cos(a) * sinh(b)
        double hre = Math.sin(Re) * Math.cosh(Im);
        double him = Math.cos(Re) * Math.sinh(Im);
        return new complex(hre, him);
    }

    public complex cos() {
//        实部 = cos(a) * cosh(b)
//        虚部 = -sin(a) * sinh(b)
        return new complex(Math.cos(Re) * Math.cosh(Im), -Math.sin(Re) * Math.sinh(Im));
    }

    public complex tan() {
//        实部 = (sin(2a) / (cos(2a) + cosh(2b)))
//        虚部 = (sinh(2b) / (cos(2a) + cosh(2b)))
        return new complex(Math.sin(2 * Re) / (Math.cos(2 * Re) + Math.cosh(2 * Im)),
                Math.sinh(2 * Im) / (Math.cos(2 * Re) + Math.cosh(2 * Im)));
    }

    public complex sinh() {
//        实部 = sinh(a) * cos(b)
//        虚部 = cosh(a) * sin(b)
        return new complex(Math.sinh(Re) * Math.cos(Im), Math.cosh(Re) * Math.sin(Im));
    }

    public complex cosh() {
//        实部 = cosh(a) * cos(b)
//        虚部 = sinh(a) * sin(b)
        return new complex(Math.cosh(Re) * Math.cos(Im), Math.sinh(Re) * Math.sin(Im));
    }

    public complex tanh(){
//        实部 = sinh(2a) / (cosh(2a) + cos(2b))
//        虚部 = sin(2b) / (cosh(2a) + cos(2b))
    return  new complex(Math.sinh(2*Re)/(Math.cosh(2*Re) + Math.cos(2*Im)),
                Math.sin(2*Im)/(Math.cosh(2*Re)+Math.cos(2*Im)));
    }


    public double getRe() {
        return Re;
    }

    public double getIm() {
        return Im;
    }

    @Override
    public String toString() {
        return "( " +
                Re + " , " +
                Im +
                " )";
    }
}

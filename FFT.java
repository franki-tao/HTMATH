package 结构性存款;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import java.util.Arrays;

/**
 * 基于apache的快速傅里叶变换，即其逆变换
 */
public final class FFT {
    public Complex[] cdata;
    private final FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
    public FFT(Complex[] cdata) {
        int n = cdata.length;
        int count = 1;
        while (n>count) {
            count*=2;
        }
        this.cdata = new Complex[count];
        System.arraycopy(cdata, 0, this.cdata, 0, n);
        for (int i = n; i < count; i++) {
            this.cdata[i] = new Complex(0,0);
        }
    }

    public Complex[] fft(){
        return fft.transform(cdata, TransformType.FORWARD);
    }

    public Complex[] ifft(){
        return fft.transform(cdata, TransformType.INVERSE);
    }

    public static void main(String[] args) {
        Complex[] cmp = new Complex[5];
        cmp[0] = new Complex(1,2);
        cmp[1] = new Complex(3,4);
        cmp[2] = new Complex(5,6);
        cmp[3] = new Complex(7,8);
        cmp[4] = new Complex(1,1);
        FFT fft = new FFT(cmp);
        System.out.println(Arrays.toString(fft.fft()));
        System.out.println(Arrays.toString(fft.ifft()));
    }
}

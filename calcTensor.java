/*Coded by Patrick O'Banion July 2017
Code adapted from matlab code "Tensor Voting Framework" by Trevor Linton
https://www.mathworks.com/matlabcentral/fileexchange/21051-tensor-voting-framework
 */
package com.example.tensor;

import android.content.Intent;
import android.content.res.Resources;
import android.graphics.Bitmap;
import android.graphics.Color;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.widget.ImageView;
import android.widget.TextView;

import java.util.ArrayList;

public class calcTensor extends AppCompatActivity {
    private ArrayList<ArrayList<Float>> dList = new ArrayList<ArrayList<Float>>();
    String dString;
    public int sigma = 20;

    @Override
    protected void onCreate(Bundle savedInstanceState) {

        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_calctensor);

        //get data set name
        Intent intent = getIntent();
        int id = intent.getIntExtra(MainActivity.id, 0);

        dString = getResources().getResourceName(id);
        dString = dString.replace("com.example.tensor:id/","");//TODO make more dynamic, search for last /

        id = getResources().getIdentifier(dString, "string", calcTensor.this.getPackageName());

        TextView edgeSet = (TextView) findViewById(R.id.edgeSet);
        edgeSet.setText(getResources().getString(id));


        float T[][][][] = readEdge();
        Tensor T1 = convertTensorEv(T);
        float[] imRange = {0,1};//[0]=min,[1]=max
        drawGrid(T1.l0(),0, imRange);//TODO scale graph to frame size
        T = findFeatures(T1.l0(),sigma);
        Tensor T2 = convertTensorEv(T);


        float[][] new2 = sub2DArray(T2.l0(),T2.l1());
        for(int i=0;i <new2.length; i++) {//find min and max
            for (int j = 0; j < new2[0].length; j++) {
                if(new2[i][j]<imRange[0]){
                    imRange[0]=new2[i][j];
                }else if(new2[i][j]>imRange[1]){
                    imRange[1]=new2[i][j];
                }
            }
        }
        drawGrid(new2,1, imRange);
        float[][] new3 = calcOrthoExtreme(T2, 15, Math.PI/8d);
        imRange[0] = 0;
        imRange[1] = 1;
        drawGrid(new3,2, imRange);


    }
    //round to E-4 precision, same as matlab
    public float round(double in){return (float)(Math.round(in*10000d)/10000d);}
    //subtracts values of 2 equal sized 2D-arrays
    public float[][] sub2DArray(float[][] i0, float[][] i1){
        float[][] output = new float[i0.length][i0[0].length];
        for (int i=0; i<i0.length; i++){
            for (int j=0; j<i0[0].length; j++ ){
                output[i][j] = i0[i][j]-i1[i][j];
            }
        }
        return output;
    }
    //reads in input file and outputs tensor field
    public float[][][][] readEdge() {
        String[] temp;//TODO code forces all initial saliencies to 1
        //EditText dArray = (EditText)findViewById(R.id.arrayText);
        Resources res = getResources();
        double x = 0, y = 0;
        int size = 0;
        int h = 0, w = 0;
        int id = 0;
        //end objects

        //get array
        id = getResources().getIdentifier(dString, "array", calcTensor.this.getPackageName());
        String[] edgeD = res.getStringArray(id);
        size = Integer.parseInt(edgeD[0]);//get size of input data

        //dArray.append("Size = " + Integer.toString(size)+"\n");//append size
        for (int i = 0; i < size; i++) {
            temp = (edgeD[i + 1]).split(" ");
            ArrayList<Float> row = new ArrayList<Float>();
            for (int g = 0; g < 3; g++) {
                row.add(Float.parseFloat(temp[g]));
            }
            dList.add(row);

            //dArray.append("\n"+String.valueOf(dList.get(i)));//Print out array
        }//TODO need to add error handling for if there is a blank line
        //find max of row 1 and row 2 separately
        for (int i = 0; i < size; i++) {
            if (dList.get(i).get(0) > h) {
                h = (int) (float) (dList.get(i).get(0));
            }
            if (dList.get(i).get(1) > w) {
                w = (int) (float) (dList.get(i).get(1));
            }
        }
        //make T array
        float[][][][] T = new float[h][w][2][2];
        for (int i = 0; i < size; i++) {//TODO round values to save space/processing time/power consumption
            x = (float) Math.cos(dList.get(i).get(2) * Math.PI / 180 + 90 * Math.PI / 180);
            y = (float) Math.sin(dList.get(i).get(2) * Math.PI / 180 + 90 * Math.PI / 180);
            T[(int) (h - dList.get(i).get(0))][(int)(dList.get(i).get(1) - 1)][0][0] = (float)(x * x);
            T[(int) (h - dList.get(i).get(0))][(int)(dList.get(i).get(1) - 1)][0][1] = (float)(x * y);
            T[(int) (h - dList.get(i).get(0))][(int)(dList.get(i).get(1) - 1)][1][0] = (float)(x * y);
            T[(int) (h - dList.get(i).get(0))][(int)(dList.get(i).get(1) - 1)][1][1] = (float)(y * y);
        }

        return T;
    }
    //convert from tensor field
    public Tensor convertTensorEv(float[][][][] T) {
        int h = T.length;
        int w = T[0].length;

        float[][][] o0 = new float[h][w][2];
        float[][][] o1 = new float[h][w][2];
        float[][] o2 = new float[h][w];
        float[][] o3 = new float[h][w];
        float a, b, t;

        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                t = (T[i][j][0][0] + T[i][j][1][1]) / 2;
                a = T[i][j][0][0] - t;
                b = (float) Math.sqrt((a * a) + (T[i][j][0][1] * T[i][j][0][1]));

                o2[i][j] = b + t;
                o3[i][j] = -b + t;

                t = (float) Math.atan2(b - a, T[i][j][0][1]);
                o0[i][j][0] = (float) Math.cos(t);
                o0[i][j][1] = (float) Math.sin(t);
                o1[i][j][0] = (float) -Math.sin(t);
                o1[i][j][1] = (float) Math.cos(t);
            }
        }
        /*EditText dArray = (EditText) findViewById(R.id.arrayText);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                dArray.append(Float.toString(o0[i][j][0]) + "\n");
            }
        }*/
    return new Tensor(o0, o1, o2, o3);
    }
    //convert to tensor field
    public float[][][][] convertTensorEv(float[][][] i0, float[][][] i1, float[][] i2, float[][] i3){
        int h = i2.length;
        int w = i2[0].length;
        float[][][][]T = new float[h][w][2][2];
        for(int i=0; i<h; i++){
            for(int j=0; j<w; j++){
                T[i][j][0][0] = (float)(i2[i][j]*Math.pow(i0[i][j][0],2)+ i3[i][j]*Math.pow(i1[i][j][0],2));
                T[i][j][0][1] = (i2[i][j]*i0[i][j][0]*i0[i][j][1]+ i3[i][j]*i2[i][j]*i1[i][j][1]);
                T[i][j][1][0] = T[i][j][0][1];
                T[i][j][1][1] = (float)(i2[i][j]*Math.pow(i0[i][j][1],2) + Math.pow(i3[i][j]*i1[i][j][1],2));
            }
        }
        return T;
    }
    //draw the images
    public void drawGrid(float[][] inIm, int num, float[] range){
        int buf=0;
        Bitmap image = Bitmap.createBitmap(inIm[0].length, inIm.length, Bitmap.Config.ARGB_8888);///////d/d/d/d/d/d/
        String name = "im" + Integer.toString(num);
        int id = getResources().getIdentifier(name, "id", calcTensor.this.getPackageName());
        ImageView im = (ImageView)findViewById(id);

        for (int y = 0; y < inIm.length; y++) {
            for (int x = 0; x < inIm[0].length; x++) {

                buf = (int)((inIm[y][x]-range[0])/(range[1]-range[0])*255);
                buf = (255<<24) | (buf<<16) | (buf<<8) | (buf);//TODO make not fuzzy
                image.setPixel(x, y, buf);
            }
        }
        im.setImageBitmap(image);
    }
    //tensor magic
    public float[][][][] findFeatures(float[][] l0, int sigma){
        float[][][][] T;
        float[][][][][] cachedVf;
        float[][][][] sparseTf;
        float[][][][] refinedTf;
        Tensor refined;
//TODO CACHEDVF
        //cachedVf = createCachedVf(sigma);//calculate voting field at various angles

        float maxL0=0;//get norm of l0
        for(int i=0; i<l0.length; i++) {
            for (int j = 0; j < l0[0].length; j++) {
                if (maxL0 < l0[i][j]) {
                    maxL0 = l0[i][j];
                }
            }
        }for(int i=0; i<l0.length; i++) {
            for (int j = 0; j < l0[0].length; j++) {
                l0[i][j] /= maxL0;
            }
        }

        //produce image of sparse tensor tokens
        sparseTf = calcSparseField(l0);
        //use ball votes
        refinedTf = calcRefinedField(sparseTf, l0, sigma);
        //clear l1
        refined = convertTensorEv(refinedTf);
        refinedTf = convertTensorEv(refined.e0(),refined.e1(),refined.l0(),new float[l0.length][l0[0].length]);

        T = calcVoteStick(refinedTf, sigma);//, cachedVf);
        return T;
    }
    //create cached Vector field from multiple angles
    public float[][][][][] createCachedVf(int sigma){
        double[] v = new double[2];//x==v[0],y==v[1]
        int wS = (int)(Math.floor(Math.ceil(Math.sqrt(-Math.log(0.01)*Math.pow(sigma,2))*2)/2)*2+1);
        float[][][][][] out = new float[180][wS][wS][2][2];
        float[][][][] fK;

        for (int i=0; i<180; i++){
            v[0] = (Math.cos(Math.PI/180.0*(i+1)));
            v[1] = (Math.sin(Math.PI/180.0*(i+1)));

            fK = createStickTf(v, sigma, wS);
            for(int a=0; a<wS; a++) {
                for (int b = 0; b < wS; b++) {
                    for (int c = 0; c < 2; c++) {
                        for (int d = 0; d < 2; d++) {
                            out[i][a][b][c][d] = fK[a][b][c][d];
                        }
                    }
                }
            }
        }
        return out;
    }

    public float[][][][] createStickTf(double[] v, int sigma, int wS){
        float[][][][]T = new float[wS][wS][2][2];//value to be returned
        double[][] rV = {{v[0],-v[1]},{v[1],v[0]}};//turn v into rotation matrix
        float bTheta = (float)Math.atan2(v[1],v[0]);
        int wHalf = (wS-1)/2;   //calculate window size from sigma

        double T0,T1;
        double s,k,l,c = (float)5.3853*(sigma-1);
        /*
        s:arc length
        k:curvature
        l:distance between voter and receiver
        c: function of the scale,
        large scale favors large interactions and makes more smooth, aid noise removal
        small scale makes more detailed
        */
        double theta;
        //angle between tangent of osculating circle at the voter and line going through voter and receiver
        double DF = 0;
        double[][] x = new double[wS][wS];
        double[][] y = new double[wS][wS];
        double[][] z = new double[2][wS*wS];

        for(int i= 0; i<wS; i++) {//TODO move rotational vector multiplication from stick to ball tensor
            for (int j = 0; j < wS; j++) {
                x[j][i] = i - wHalf;
                y[i][j] = wHalf - i;
            }
        }
        for(int i=0; i<wS; i++) {//matlab meshgrid function
            for(int j=0; j<wS; j++) {
                z[0][(i*wS)+j] = rV[0][0] * x[j][i]+rV[1][0] * y[j][i];
                z[1][(i*wS)+j] = rV[0][1] * x[j][i]+rV[1][1] * y[j][i];
            }
        }

        for(int i=0; i<wS; i++) {//set values in the T array
            for(int j=0; j<wS; j++) {
                theta = Math.atan2(z[1][(i*wS)+j],z[0][(i*wS)+j]);
                T0=-Math.sin(2d*theta+bTheta);
                T1=Math.cos(2d*theta+bTheta);

                theta = Math.abs(theta);
                if(theta>(Math.PI/2d)){
                    theta=Math.PI-theta;
                }
                theta *= 4d;

                l = Math.sqrt(x[j][i]*x[j][i]+y[j][i]*y[j][i]);
                if(l != 0 && theta != 0) {
                    s = (theta * l / Math.sin(theta));
                }else {
                    s = l;
                }if (l != 0){
                    k = (2d*Math.sin(theta)/l);
                }else
                    k = 0;
                DF = Math.exp(-((s*s+c*k*k)/(sigma*sigma)));
                if(theta>Math.PI/2d) {
                    DF = 0;
                }
                //round the answers because 0.0000004 might as well be 0
                T[j][i][0][0]= round(T0*T0*DF);
                T[j][i][0][1]= round(T0*T1*DF);
                T[j][i][1][0]= round(T1*T0*DF);
                T[j][i][1][1]= round(T1*T1*DF);

            }
        }

        return T;
    }

    public float[][][][] calcSparseField(float[][] l0){
        float sparseTf[][][][] = new float[l0.length][l0[0].length][2][2];
            for (int i=0; i<l0.length; i++){
                for (int j=0; j<l0[0].length; j++){
                    if(l0[i][j]>0){
                        sparseTf[i][j][1][1] = sparseTf[i][j][0][0] = 1;
                    }
                }
            }

        return sparseTf;
    }

    //refine voteBall to proper places
    public float[][][][] calcRefinedField(float[][][][] tF, float[][] l0, int sigma){
        float refTf[][][][] = calcVoteBall(tF, l0, sigma);
        for(int i=0; i<refTf.length; i++){
            for(int j=0; j<refTf[0].length; j++){
                if(l0[i][j]!=0){//only add values where l0 is not 0
                    tF[i][j][0][0]+=refTf[i][j][0][0];
                    tF[i][j][0][1]+=refTf[i][j][0][1];
                    tF[i][j][1][0]+=refTf[i][j][1][0];
                    tF[i][j][1][1]+=refTf[i][j][1][1];
                }
            }

        }

        return tF;
    }

    public float[][][][] calcVoteBall(float[][][][] sparseTf, float[][] l0, int sigma){
        int wS = (int)(Math.floor(Math.ceil(Math.sqrt(-Math.log(0.01)*Math.pow(sigma,2))*2d)/2d)*2d+1d);
        int wHalf = (wS-1)/2;
        int tH = l0.length;
        int tW = l0[0].length;
        float [][][][] ballTf = createBallTf(sigma, wS);
        float voteBall[][][][] = new float[tH+wHalf*2][tW+wHalf*2][2][2];
        float output[][][][] = new float[tH][tW][2][2];
        ArrayList<ArrayList<Integer>> D = new ArrayList<>();
        D.add(new ArrayList<Integer>());
        D.add(new ArrayList<Integer>());
        int x,y;

        for (int i=0; i<tH; i++){
            for( int j=0; j<tW; j++){
                voteBall[i+wHalf][j+wHalf][0][0] = sparseTf[i][j][0][0];
                voteBall[i+wHalf][j+wHalf][0][1] = sparseTf[i][j][0][1];
                voteBall[i+wHalf][j+wHalf][1][0] = sparseTf[i][j][1][0];
                voteBall[i+wHalf][j+wHalf][1][1] = sparseTf[i][j][1][1];
            }
        }
        Tensor ball = convertTensorEv(voteBall);
        for (int q=0; q<ball.l0YLength(); q++){
            for (int w=0; w<ball.l0XLength(); w++){
                if(ball.l0(q,w)>0) {
                    D.get(0).add(q);
                    D.get(1).add(w);
                }
            }
        }
        for(int p=0; p<D.get(0).size(); p++){
            x=D.get(0).get(p) - wHalf;
            y=D.get(1).get(p) - wHalf;
            for (int bX=0; bX < wHalf*2; bX++) {
                for (int  bY=0; bY<wHalf*2; bY++) {
                    voteBall[bX+x][bY+y][0][0] += ballTf[bX][bY][0][0] * l0[x][y];
                    voteBall[bX+x][bY+y][0][1] += ballTf[bX][bY][0][1] * l0[x][y];
                    voteBall[bX+x][bY+y][1][0] += ballTf[bX][bY][1][0] * l0[x][y];
                    voteBall[bX+x][bY+y][1][1] += ballTf[bX][bY][1][1] * l0[x][y];
                }
            }
        }
        for(int l=0; l<tH; l++){
            for (int w=0; w<tW; w++){
                output[l][w][0][0] = voteBall[wHalf+l][wHalf+w][0][0];
                output[l][w][0][1] = voteBall[wHalf+l][wHalf+w][0][1];
                output[l][w][1][0] = voteBall[wHalf+l][wHalf+w][1][0];
                output[l][w][1][1] = voteBall[wHalf+l][wHalf+w][1][1];
            }
        }
        return output;
    }

    public float[][][][] createBallTf(int sigma, int wS){
        float[][][][] ballTf = new float[wS][wS][2][2];
        float[][][][] buf;
        double[] v =  new double[2];
        //integration, to change how many times done, change 32 and 16
        for( int i=0; i<32; i++){
            v[0] = Math.cos(i/16d*Math.PI);
            v[1] = Math.sin(i/16d*Math.PI);
            buf = createStickTf(v, sigma, wS);
            for(int a=0; a<wS; a++){
                for(int b=0; b<wS; b++){
                    for(int c=0; c<2; c++){
                        for(int d=0; d<2; d++){
                            ballTf[a][b][c][d] += buf[a][b][c][d];
                        }
                    }
                }
            }
        }
        for(int q=0; q<wS; q++){
            for(int w=0; w<wS; w++){
                for(int e=0; e<2; e++){
                    for(int r=0; r<2; r++){
                        ballTf[q][w][e][r] = round(ballTf[q][w][e][r]/32d);
                    }
                }
            }
        }
        /*alt code to approximate integration
        remove rotational vector multiplication and instead replace with eq 2.6
        ballTf=sum from 1->k(stick tensors)(v(stick vote in vector form)*v^T)
         The stick votes from O to P cast by K stick tensors at angle intervals of 2π/K span
         the unit circle. Normalization has to be performed in order to make the energy emitted
         by a unit ball equal to that ofa unit stick
         p.26
        */
        return ballTf;
    }

    public float[][][][] calcVoteStick(float[][][][] refinedTf, int sigma){//, float[][][][][] cachedVf){

        int wS = (int)(Math.floor(Math.ceil(Math.sqrt(-Math.log(0.01)*Math.pow(sigma,2))*2d)/2d)*2d+1d);
        int wHalf = (wS-1)/2;//TODO clean up
        int tH = refinedTf.length;
        int tW = refinedTf[0].length;
        float voteStick[][][][] = new float[tH+wHalf*2][tW+wHalf*2][2][2];
        ArrayList<ArrayList<Integer>> D = new ArrayList<>();
        D.add(new ArrayList<Integer>());
        D.add(new ArrayList<Integer>());
        int x,y;
        double[] v = new double[2];
        float [][][][] stickTf;
        float output[][][][] = new float[tH][tW][2][2];
        float weight;

        for (int i=0; i<tH; i++){
            for( int j=0; j<tW; j++){
                voteStick[i+wHalf][j+wHalf][0][0] = refinedTf[i][j][0][0];
                voteStick[i+wHalf][j+wHalf][0][1] = refinedTf[i][j][0][1];
                voteStick[i+wHalf][j+wHalf][1][0] = refinedTf[i][j][1][0];
                voteStick[i+wHalf][j+wHalf][1][1] = refinedTf[i][j][1][1];
            }
        }
        Tensor stick = convertTensorEv(voteStick);

        for (int q=0; q<stick.l0YLength(); q++){
            for (int w=0; w<stick.l0XLength(); w++){
                if(stick.l0(q,w)-stick.l1(q,w)>0) {
                    D.get(0).add(q);
                    D.get(1).add(w);
                }
            }
        }
        for(int p=0; p<D.get(0).size(); p++){
            v[0]= -stick.e0(D.get(0).get(p),D.get(1).get(p),1);
            v[1]= stick.e0(D.get(0).get(p),D.get(1).get(p),0);
            stickTf = createStickTf(v, sigma, wS);
            //get weight for stickTf
            weight = stick.l0(D.get(0).get(p),D.get(1).get(p))-stick.l1(D.get(0).get(p),D.get(1).get(p));

            x=D.get(1).get(p) - wHalf;
            y=D.get(0).get(p) - wHalf;
            for (int bX=0; bX < wHalf*2; bX++) {
                for (int  bY=0; bY<wHalf*2; bY++) {
                    voteStick[bY+y][bX+x][0][0] += (stickTf[bY][bX][0][0] * weight);
                    voteStick[bY+y][bX+x][0][1] += (stickTf[bY][bX][0][1] * weight);
                    voteStick[bY+y][bX+x][1][0] += (stickTf[bY][bX][1][0] * weight);
                    voteStick[bY+y][bX+x][1][1] += (stickTf[bY][bX][1][1] * weight);
                }
            }
        }

        for(int l=0; l<tH; l++){
            for (int w=0; w<tW; w++){
                output[l][w][0][0] = voteStick[wHalf+l][wHalf+w][0][0];
                output[l][w][0][1] = voteStick[wHalf+l][wHalf+w][0][1];
                output[l][w][1][0] = voteStick[wHalf+l][wHalf+w][1][0];
                output[l][w][1][1] = voteStick[wHalf+l][wHalf+w][1][1];
            }
        }
        return output;
    }

    //finds orthogonal of the curve and maxima along the line
    // r is length to check along the normal, epsilon is width in radians to check
    public float[][] calcOrthoExtreme(Tensor tens, int r, double epsilon){
        int buf = r*2+1;
        double bufd;
        float[][] q = sub2DArray(tens.l0(),tens.l1());
        float[][] q1 = new float[2*r+tens.l0YLength()][2*r+tens.l0XLength()];
        float[][] out = new float[tens.l0YLength()][tens.l0XLength()];
        float[][] x = new float[buf][buf];
        float[][] y = new float[buf][buf];
        float[] z = new float[buf];
        float[][] t = new float[buf][buf];
        float[][] l = new float[buf][buf];
        //double x2, y2, t2;
        double[][] x2 = new double[buf][buf];
        double[][] y2 = new double[buf][buf];
        double[][] t2 = new double[buf][buf];
        boolean val;
        ArrayList<ArrayList<Integer>> D = new ArrayList<>();
        D.add(new ArrayList<Integer>());
        D.add(new ArrayList<Integer>());

        for (int i=(0-r); i<r+1; i++) {
            for(int j=(0-r); j<r+1; j++) {
                x[i+r][j+r] = j;
                y[j+r][i+r] = j;
            }
        }
        for (int g=0; g<x.length; g++){
            for (int w=0; w<x[0].length; w++){
                t[g][w] = (float)Math.atan2(y[g][w],x[g][w]);
                l[g][w] = (float)Math.sqrt(x[g][w]*x[g][w]+y[g][w]*y[g][w]);
            }
        }
        for (int a=r; a<tens.l0YLength()+r; a++){
            for (int b=r; b<tens.l0XLength()+r; b++){
                q1[a][b] = q[a-r][b-r];
                if(q1[a][b] > 0){
                    D.get(0).add(a);//y
                    D.get(1).add(b);//x
                }
            }
        }

        for(int e=0; e<D.get(0).size(); e++){
            val = true;
            bufd = Math.atan2(tens.e0(D.get(0).get(e)-r,D.get(1).get(e)-r,1),tens.e0(D.get(0).get(e)-r,D.get(1).get(e)-r,0));
            for(int f=0; f<l.length; f++){//y
                for(int g=0; g<l[0].length; g++) {//x
                    //if(q1[D.get(0).get(e)][D.get(1).get(e)] != 0 && l[f][g]<=r) {
                    if(l[f][g] <= r){
                        x2[f][g] = l[f][g] * Math.cos(t[f][g] + bufd);
                        y2[f][g] = l[f][g] * Math.sin(t[f][g] + bufd);
                        t2[f][g] = Math.abs(Math.atan2(y2[f][g], x2[f][g]));
                        if (t2[f][g] > Math.PI / 2d) {
                            t2[f][g] = Math.PI - t2[f][g];
                        }
                        if (t2[f][g] <= epsilon) {
                            t2[f][g] = 1;
                            /*if(q1[f][g]>q1[D.get(0).get(e)][D.get(1).get(0)]){
                                val = false;
                            }
                            //if (q1[f + r][g + r] < q1[D.get(0).get(e)][D.get(1).get(e)]) {
                              //  val = false;
                            //}*/
                        }if(t2[f][g] != 1){
                            t2[f][g] = 0;
                        }
                    }
                }
            }for(int c=0; c<l.length; c++){
                for(int d=0; d<l[0].length; d++){
                    if(q1[D.get(0).get(e)-r+c][D.get(1).get(e)-r+d]*t2[c][d]>q1[D.get(0).get(e)][D.get(1).get(e)]){
                        val = false;
                    }
                }

            }
            if(val){
                out[D.get(0).get(e)-r][D.get(1).get(e)-r] = 1;
            }
        }
        return out;
    }
}

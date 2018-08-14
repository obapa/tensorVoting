package com.example.tensor;

/**
 * Basic class defining Tensor object, used in calcTensor
 */


public class Tensor{
    private float[][][] e0,e1;
    private float[][] l0,l1;

    public Tensor(float[][][] e0, float[][][] e1, float[][] l0, float[][] l1){
        this.e0=e0;
        this.e1=e1;
        this.l0=l0;
        this.l1=l1;
    }


    //gets
    public float[][][] e0(){ return e0;}
    public float[][][] e1(){ return e1;}
    public float[][] l0(){ return l0;}
    public float[][] l1(){ return l1;}
    public float e0(int i, int j, int g){ return e0[i][j][g];}
    public float e1(int i, int j, int g){ return e1[i][j][g];}
    public float l0(int i, int j){ return l0[i][j];}
    public float l1(int i, int j){ return l1[i][j];}
    public int l0YLength() {return l0.length;}
    public int l0XLength() {return l0[0].length;}


    //sets
    public void sete0(float[][][]e0){this.e0=e0;}
    public void sete1(float[][][]e1){this.e1=e1;}
    public void setl0(float[][]l0){this.l0=l0;}
    public void setl1(float[][]l1){this.l1=l1;}

}

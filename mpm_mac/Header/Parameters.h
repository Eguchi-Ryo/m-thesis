#pragma once
//プログラム全般
//static double DT = 1.0/ 12000;
const static double BSPLINE_EPSILON = 1e-8;
const static double MASS_EPSILON = 1e-8;
const static double CFRI = 0.2;
const static double DAMPRATIO = 0.8;


//描画用
const static int WINDOW_WIDTH = 700; //ウィンドウの幅
const static int WINDOW_HEIGHT = 700; //ウィンドウの高さ
const static double VIEW_WIDTH = 1.0;
const static double VIEW_HEIGHT = WINDOW_HEIGHT * VIEW_WIDTH / WINDOW_WIDTH;




#define ORGAN_KIND 4
#define INIT_THRESHOLD 0.005
//十二指腸の弾性体パラメータ
const static double DUO_MIU = 3.38e6;
const static double DUO_LAMBDA = 8.11e7;



//膵臓の弾性体パラメータ
const static double PANC_MIU = 1.17e3;
const static double PANC_LAMBDA = 5.76e4;

//腎臓の弾性体パラメータ
const static double LKID_MIU = 8.01e3;
const static double LKID_LAMBDA = 3.99e6;

//胃の弾性体パラメータ
// young 500kPa(何か人が食べていることを考慮。しない場合は1~100kPa程度)
// poisson 0.499
const static double STOM_MIU = 1.68e5;
const static double STOM_LAMBDA = 8.32e7;










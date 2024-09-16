#include <iostream>
#include <vector>
#include <random>
#include <cmath>

using namespace std;

int main()
{
    // Constant Definition
    int L = 100;        // L*L 격자 실험
    int epochs = 1000;  // 시뮬레이션 횟수
    double J = 1.0;
    double kT = 2.2;       // 실험 설정 온도 / [1] critical temperature (Tc)보다 약간 낮게 설정한다.

    //Random number Generation
    random_device rd;         // 시드값을 얻기 위한 random_device 생성.
    mt19937 gen(rd());        // random_device(rd()) 를 통해 난수 생성 엔진을 초기화 - mt19937 Mersenne Twister random number Generator
    uniform_int_distribution<int> dis(0, 1);    // 0과 1 까지 균등하게 나타나는 난수열을 생성하기 위해 균등 분포 정의.

    //Initialize State
    // -1 or 1의 랜덤 숫자로 구성된 L*L 격자 구조 State 생성
    vector<vector<int>> state;

    for(int i=0;i<L;i++)
    {
        vector<int> num;
        int ran;
        for(int j=0;j<L;j++)
        {
            ran = dis(gen);
            if(ran) num.push_back(ran);
            else    num.push_back(-1);
        }
        state.push_back(num);
        num.clear();
    }

    // Ising Model Monte Carlo Method - Metropolis Algorithms

    int x, y;   // 임의의 State 위치 변수
    int state_xy, nbstate;  // state 값들 저장 변수
    int dH;     // 변경 전후, 에너지 차이값

    double prob;
    uniform_real_distribution<double> distr(0,1);

    for(int rep=0;rep<epochs;rep++)
    {
        // 임의의 State 선택
            //uniform_int_distribution<int> dis_x(0, L);
        x = dis(gen)*L;
        y = dis(gen)*L;
        
        state_xy = state[x][y];    //임의의 State 선택
        
        // nearest-neighbor states 합계 계산
        nbstate = state[(x+1)%L][y] + state[(x-1)%L][y] + state[x][(y+1)%L] + state[x][(y-1)%L];
        
        // Flip 여부 판단
        /*
            S의 state를 flip해 보고 dE에 따라 flip할 것인지 결정한다 [2].
            dE < 0이면 flip을 허용하고, dE > 0이면 exp(-dE / kT)의 확률로 flip을 허용한다.
                H1 = -J * nbState * S
                H2 = -J * nbState * (S * (-1))
                dH = H2 - H1
        */

        dH = J * nbstate * state_xy * 2;

        if(dH < 0)
            state[x][y] *=-1;            //flip spin
        else
        {
            prob = exp(-(double)dH / kT);
            if (distr(gen)<= prob)
                state[x][y] *= -1;       // flip spin
        }
    }

    // Calculate Energy and Megnetization


    // Storing Data


    return 0;
}
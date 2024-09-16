#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <random>
#include <cmath>

using namespace std;

int main()
{
    // Constant Definition
    int L = 30;        // L*L 격자 실험
    int eq_steps = 10000;   // equilibrium value가 될때까지의 시뮬레이션 반복 횟수
    int mc_steps = 1000;    // equilibrium value에 도달한 후, 시뮬레이션 결과값 계산을 위한 횟수 (평균, 분산 등)
    double J = 1.0;
    double interval = 0.01; // 측정 온도의 간격
    double end_kT = 4.0;
    //double kT = 2.2;       // 실험 설정 온도 / [1] critical temperature (Tc)보다 약간 낮게 설정한다.

    vector<double> kT_range;
    for(double k = interval; k< end_kT;k+=interval)   kT_range.push_back(k);

    //Random number Generation
    random_device rd;         // 시드값을 얻기 위한 random_device 생성.
    mt19937 gen(rd());        // random_device(rd()) 를 통해 난수 생성 엔진을 초기화 - mt19937 Mersenne Twister random number Generator
    uniform_int_distribution<int> dis(0, 1);    // 0과 1 까지 균등하게 나타나는 난수열을 생성하기 위해 균등 분포 정의.
    uniform_int_distribution<int> dis_l(0, L-1);    // 0과 1 까지 균등하게 나타나는 난수열을 생성하기 위해 균등 분포 정의.

    // range kT 온도별 시뮬레이션 반복
    vector<double> E, M, C, X;

    cout << L<<"x"<<L <<"Start\n";
    for(double kT : kT_range)
    {
        cout << kT/interval <<"/"<< end_kT/interval<<"\n";

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
        
        for(int rep=0;rep<eq_steps;rep++)
        {
            // 임의의 State 선택
                //uniform_int_distribution<int> dis_x(0, L);
            x = dis_l(gen);
            y = dis_l(gen);
            
            state_xy = state[x][y];    //임의의 State 선택

            // nearest-neighbor states 합계 계산
            if(x==0 & y==0)
                nbstate = state[(x+1)%L][y] + state[L-1][y] + state[x][(y+1)%L] + state[x][L-1];
            else if(x==0)
                nbstate = state[(x+1)%L][y] + state[L-1][y] + state[x][(y+1)%L] + state[x][y-1];
            else if(y==0)
                nbstate = state[(x+1)%L][y] + state[x-1][y] + state[x][(y+1)%L] + state[x][L-1];
            else
                nbstate = state[(x+1)%L][y] + state[x-1][y] + state[x][(y+1)%L] + state[x][y-1];
            
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

        // MC_step energy / magnetization 값 계산을 위한 평균/분산 값
        double E1= 0, M1=0, E2 = 0, M2 = 0;
        for(int rep=0;rep<mc_steps;rep++)
        {
          // 임의의 State 선택
                //uniform_int_distribution<int> dis_x(0, L);
            x = dis_l(gen);
            y = dis_l(gen);
            
            state_xy = state[x][y];    //임의의 State 선택

            // nearest-neighbor states 합계 계산
            if(x==0 & y==0)
                nbstate = state[(x+1)%L][y] + state[L-1][y] + state[x][(y+1)%L] + state[x][L-1];
            else if(x==0)
                nbstate = state[(x+1)%L][y] + state[L-1][y] + state[x][(y+1)%L] + state[x][y-1];
            else if(y==0)
                nbstate = state[(x+1)%L][y] + state[x-1][y] + state[x][(y+1)%L] + state[x][L-1];
            else
                nbstate = state[(x+1)%L][y] + state[x-1][y] + state[x][(y+1)%L] + state[x][y-1];
            
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
  


            // Calculate Energy and Megnetization
            int energy = 0;
            int mag = 0;
            for(int i=0;i<L;i++)
            {
                for(int j=0;j<L;j++)
                {
                    // nearest-neighbor states 합계 계산
                    if(x==0 & y==0)
                        nbstate = state[(x+1)%L][y] + state[L-1][y] + state[x][(y+1)%L] + state[x][L-1];
                    else if(x==0)
                        nbstate = state[(x+1)%L][y] + state[L-1][y] + state[x][(y+1)%L] + state[x][y-1];
                    else if(y==0)
                        nbstate = state[(x+1)%L][y] + state[x-1][y] + state[x][(y+1)%L] + state[x][L-1];
                    else
                        nbstate = state[(x+1)%L][y] + state[x-1][y] + state[x][(y+1)%L] + state[x][y-1];

                    energy += -nbstate*state[i][j];
                    mag += state[i][j];
                }
            }
            energy /= 4;    // nn 합계 후, 중복 계산된 것들을 나눔

            E1 += energy;
            M1 += mag;
            E2 += energy*energy;
            M2 += mag*mag;
        }

        double n1 = 1.0/(mc_steps*L*L);
        double n2 = 1.0/(mc_steps*mc_steps*L*L);
        double iT = 1.0/kT, iT2 = 1.0/(kT*kT);

        E.push_back(E1*n1);
        M.push_back(M1/n1);
        C.push_back((n1*E2 - n2*E1*E1)*iT2);
        X.push_back((n1*M2 - n2*M1*M1)*iT);

        //state 초기화
        state.clear();
    }

    // Storing Data
    string filePATH = "./dataset/data_"+ to_string(L) +"x"+to_string(L)+"_eqsteps"+to_string(eq_steps)+".txt";

    ofstream writeFile(filePATH.data());
    if(writeFile.is_open())
    {
        writeFile << "Temperature(kT) \t Energy \t Magnetization \t Specific Heat \t Susceptibility \n";

        for(int i =0; i< kT_range.size(); i++)
            writeFile << to_string(kT_range[i]) << "\t" << to_string(E[i]) << "\t" << to_string(M[i]) << "\t" << to_string(C[i]) << "\t" << to_string(X[i]) << "\n";

        writeFile.close();
    }
    cout << "End";


    return 0;
}
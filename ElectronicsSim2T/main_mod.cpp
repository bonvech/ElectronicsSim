/**
 * Версия 1.0.2 от 13.12.2023
 */
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <sstream>
#include <numeric>
#include <thread>
#include <functional>
#include <random>
#include <cstring>


#define PULSE_LENGTH 1000
#define N_CHAN 109    /// Число каналов    
#define BIN_2_GEN 512
#define AMP_SIZE 10000   /// Длина импульса электртоники
#define INTERF_LENGTH 6200 /// Длина наводки


static float CURR_2_PH = 3. / 8;
static float pieds[2*N_CHAN];      /// Пьедесталы
static float curbase[N_CHAN];    /// Относительные токи
static float pulse[PULSE_LENGTH];      /// Импульсные характеристики тока
static float amp[AMP_SIZE];        /// Обратная функция распределения коэффициента усиления ФЭУ
static float Toff[N_CHAN];       /// Относительные сдвиги каналов
static float interf [INTERF_LENGTH];     /// Наводки
static float interf_amp[N_CHAN]; /// Покональные амплитудные коэффициенты наводки
static int thr[N_CHAN];          /// Пороги
static int triggermask[N_CHAN];  /// Маска триггера
static float Cal[N_CHAN];        /// Калибровка
static int N_PHEL = 0;                  /// Число зарегистрированных фотонов
std::vector<int> PMTid;        /// Номера ФЭУ, куда упали фотоэлектроны
std::vector<float> T;          /// Времена прихода фотоэлектронов
static float Tmin{};                  /// Минимальное время прихода
static float MEAN_CURR = 3.5;         /// Средний ток
static int data_out[BIN_2_GEN+1][N_CHAN]; /// Выводной массив
static float data [N_CHAN][BIN_2_GEN * 50 + 2 * PULSE_LENGTH + 2]; /// Данные

static int BG_LENGTH{50 * BIN_2_GEN + PULSE_LENGTH + 50};
static float N_AVG, N_PHEL_exp, amp_ph;
static int T_ph;

static int L2_links[N_CHAN][6]={0};  // в L2 у каждого ФЭУ только 6 связей
static int L2_logic[290]={0};

bool TRIGGER={0};

int i,j,TCOUNT;
static int STAT = 100;



int main() {

//    manager.inputAll();	
	std::string line;
	
	
//   ЗАГРУЗКА ДАННЫЗ ФОТОЭЛЕКТРОНОВ ИЗ moshits	
	{std::ifstream moshits("mosaic_hits");
    if (!moshits.is_open()) {
        std::cerr << "Failed to open the moshits file!" << std::endl;
    }
    N_PHEL = int(std::count(std::istreambuf_iterator<char>(moshits),
                            std::istreambuf_iterator<char>(), '\n')) - 1;

    moshits.clear();
    moshits.seekg(0, std::ios::beg);

	PMTid.reserve(N_PHEL);
	T.reserve(N_PHEL);

    std::getline(moshits, line);
    double pmt_tmp, t_tmp, tmp;
    while (std::getline(moshits, line)) {
        std::istringstream ss(line);
        ss >> pmt_tmp;
        for (int i{0}; i < 3; i++) {
            ss >> tmp;
        }
        ss >> t_tmp;
        PMTid.push_back(int(pmt_tmp));
        T.push_back(float(t_tmp));
    }
    auto min_it{std::min_element(T.begin(), T.end())};
    Tmin = *min_it;
    moshits.close();
	}

//    ЗАГРУЗКА КОЭФФИЦИЕНТОВ КАЛИБРОВКИ
	float tmpf;
	int tmpi;
	{std::ifstream cal("Current.cal");
		if (!cal.is_open()) {
			std::cerr << "Failed to open the calibration file!" << std::endl;
		}		
		j = 0;
		while (std::getline(cal, line) && j < N_CHAN ) {
			std::istringstream ss(line);
			for (int i{0}; i < 13; i++) {
				ss >> tmpf;
			}
			Cal[j] = float (tmpf);
			j++;
		}
		cal.close();
	}
	
//   ЗАГРУЗКА ОТНОСИТЕЛЬНЫХ КОЭФФИЦИЕНТОВ ФОНОВЫХ ТОКОВ
	{std::ifstream input("CurRels.dat");
		if (!input.is_open()) {
			std::cerr << "Failed to open the CurRels.dat file!" << std::endl;
		}		
		std::getline(input, line);
		std::istringstream ss(line);    		
		i = 0;
		while (ss >> tmpf && i < N_CHAN) {
			curbase[i]= float (tmpf);
			i++;
		}
		input.close();
	}
	
//    ЗАГРУЗКА ФУНКЦИИ РАСПРЕДЕЛЕНИЯ КОЭФФИЦИЕНТА УСИЛЕНИЯ ФЭУ
	{std::ifstream input("AmpDistrib.dat");
		if (!input.is_open()) {
			std::cerr << "Failed to open the AmpDistrib.dat file!" << std::endl;
		}		
		std::getline(input, line);
		std::istringstream ss(line);    		
		i = 0;
		while (ss >> tmpf && i < AMP_SIZE) {
			amp[i]= float (tmpf);
			i++;
		}
		input.close();
	}
	
//     ЗАГРУЗКА ПЬЕДЕСТАЛОВ	
	{std::ifstream input("pieds.dat");
		if (!input.is_open()) {
			std::cerr << "Failed to open the pieds.dat file!" << std::endl;
		}		
		std::getline(input, line);
		std::istringstream ss(line);    		
		i = 0;
		while (ss >> tmpf && i < 2*N_CHAN) {
			pieds[i]= float (tmpf);
			i++;
		}
		input.close();
	}
	
//   ЗАГРУЗКА ОТНОСИТЕЛЬНЫХ СДВИГОВ КАНАЛОВ	
	{std::ifstream input("Toff.dat");
		if (!input.is_open()) {
			std::cerr << "Failed to open the Toff.dat file!" << std::endl;
		}		
		std::getline(input, line);
		std::istringstream ss(line);    		
		i = 0;
		while (ss >> tmpf && i < N_CHAN) {
			Toff[i]= float (tmpf);
			i++;
		}
		input.close();
	}
	
//     ЗАГРУЗКА ПОРОГОВ	
	{std::ifstream input("thresholds.dat");
		if (!input.is_open()) {
			std::cerr << "Failed to open the thresholds.dat file!" << std::endl;
		}		
		std::getline(input, line);
		std::istringstream ss(line);    		
		i = 0;
		while (ss >> tmpi && i < N_CHAN) {
			thr[i]= int (tmpi);
			i++;
		}
		input.close();
	}
	
//     ЗАГРУЗКА МАСКИ ТРИГГЕРА
	{std::ifstream input("trigger.dat");
		if (!input.is_open()) {
			std::cerr << "Failed to open the trigger.dat file!" << std::endl;
		}		
		std::getline(input, line);
		std::istringstream ss(line);    		
		i = 0;
		while (ss >> tmpi && i < N_CHAN) {
			triggermask[i]= int (tmpi);
			i++;
		}
		input.close();
	}
	
//     ЗАГРУЗКА ЛОГИКИ ТРИГГЕРА L2
	{std::ifstream input("SPHERE2_trigger_L2.dat");
		if (!input.is_open()) {
			std::cerr << "Failed to open the SPHERE2_trigger_L2.dat file!" << std::endl;
		j = 0;
		while (std::getline(input, line) && j < N_CHAN ) {
			std::istringstream ss(line);
			for (int i{0}; i < 6; i++) {
				ss >> tmpi;
				L2_links[j][i] = int (tmpi);
			}
			j++;
		}
		input.close();
		}
	}
	
//	model.GetImpulseData();
	{
		std::ifstream input("Impulse2GHz.dat");
		if (!input.is_open()) {
			std::cerr << "Failed to open the Impulse2GHz.dat file!" << std::endl;
		}
		std::string line;
		std::getline(input, line);
		std::istringstream ss(line);
		float tmpf;
		i = 0;
		while (ss >> tmpf && i< PULSE_LENGTH) {
			pulse[i]= float (tmpf);
			i++;        
		}
		input.close();
	}
	
	
	
//    НАЧАЛО ЦИКЛА ПО РЕАЛИЗАЦАЯМ

    for ( int NREAL = 0; NREAL < STAT; NREAL++){
		
// 		СБРОС МАССИВОВ ДАННЫХ НА 0
		memset(data, 0, sizeof(data));
		memset(data_out, 0, sizeof(data_out));
		TRIGGER = 0;
		
		
//    model.AddBackground();
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dist_amp(0, AMP_SIZE);
		std::uniform_int_distribution<> dist_bg(0, BG_LENGTH);

		for (int j{0}; j < N_CHAN; j++) {
			if (!triggermask[j]) {
				continue;
			}
			N_AVG = MEAN_CURR * curbase[j] * CURR_2_PH;
			N_PHEL_exp = N_AVG * float(BG_LENGTH) / 50;
			std::poisson_distribution<> dist(N_PHEL_exp);
			int N_BG_PHEL = dist(gen);
			for (int n{0}; n < N_BG_PHEL; n++) {
				amp_ph = amp[dist_amp(gen)];
				T_ph = dist_bg(gen);
				for (int t{0}; t < PULSE_LENGTH; t++) {
					data[j][t + T_ph] += amp_ph * pulse[t];
				}
			}
		}
		

//    model.SubtractMeans();
		
		float S_avg{0};
		for (int j{0}; j < N_CHAN; j++) {
			if (!triggermask[j]) {
				continue;
			}
			S_avg=0;
			for (int t{PULSE_LENGTH}; t < BIN_2_GEN * 50 + PULSE_LENGTH + 1; t++) {
				S_avg += data[j][t];
			}
			S_avg /= float(BIN_2_GEN) * 50;
			for (int t{0}; t < BIN_2_GEN * 50 + PULSE_LENGTH + 1; t++) {
				data[j][t]-=S_avg;
			}
		}
		
		
//    model.GenerateEvent();
				
		std::uniform_int_distribution<> dist(0, AMP_SIZE);
		for (int phid{0}; phid < N_PHEL; phid++) {			
			amp_ph = amp[dist(gen)];
			T_ph = int(2 * (T[phid] - Tmin) + floor( 240./512. * BIN_2_GEN * 50 + PULSE_LENGTH));
			for (int t{0}; t < PULSE_LENGTH; t++) {
				data[PMTid[phid]][t + T_ph] += amp_ph * pulse[t];
			}
		}
		
		
//    model.SimulateDig();
		
		std::uniform_int_distribution<> dist_inter_length(0, INTERF_LENGTH);
		std::uniform_int_distribution<> dist_shift(0, 50);

		int t_f{0};
		int t_shift{dist_shift(gen)};
		t_f = dist_inter_length(gen);
		t_shift = dist_shift(gen);
		for (int j{0}; j < N_CHAN; j++) { 
			if (!triggermask[j]) {
				continue;
			}		
			data_out[0][j]=int(Toff[j]);
			for (int i{0}; i < BIN_2_GEN; i++) {
				data_out[i+1][j] = int((data[j][PULSE_LENGTH + t_shift + i * 50 + int(25*Toff[j]/100.)]) / Cal[j]
					+ pieds[j*2+1]   // на триггер идут только второй субканал (причины сугубо исторические и приборные
					+ interf_amp[j] * interf[(i * 50 + int(25*Toff[j]/100.) + t_f) % INTERF_LENGTH]);
			}
		}
		
//    model.PrintDataOut();
	/*
		std::ofstream outFile("data_out");
		if (!outFile.is_open()) {
			std::cerr << "Open file error." << std::endl;
		}
		for ( i = 0 ; i < BIN_2_GEN +1 ; i++ ){
			for ( j = 0; j < N_CHAN ; j++){
			if (data_out[i][j] < 10) {outFile <<  ' ';};
			if (data_out[i][j] < 100) {outFile << ' ';};
			outFile << data_out[i][j] << ' ';
			}
			outFile << std::endl;
		}
		outFile.close();
	*/

//    model.TriggerCheck();
		
		int DiscInput[N_CHAN]={0};
		int Timers[N_CHAN]={0};
		int Tstart=50;            // начало просмотра кадра на предмет срабатывания триггера, 
		int i,j;
		
		for (j = 0 ; j < N_CHAN ; j++) {
			for (i = Tstart ; i < (4+Tstart) ; i++ ){
				DiscInput[j]+=data_out[i][j];
			}
		}
		
		i = Tstart + 4;
		while ((i < BIN_2_GEN )&&(!TRIGGER)){
			i++;
			memset(L2_logic, 0, sizeof(L2_logic));
			for ( j = 0 ; j < N_CHAN ; j++){
				if (!triggermask[j]) {
					continue;
				}				
				DiscInput[j] += (data_out[i][j] - data_out[i-4][j]);
				if (DiscInput[j] > thr[j]) {Timers[j] = 41;}
				if (Timers[j] > 0) {
					Timers[j]--;
					for (int zi = 0; zi < 6 ; zi++){
						L2_logic[L2_links[j][zi]]++;
					}
				}
			}
			if ((*std::max_element(std::begin(L2_logic), std::end(L2_logic))) > 2){		
//				std::cout << NREAL << "\t" << i << std::endl;
				TRIGGER = 1;
				TCOUNT++;
			}
		}		
	}
	
	std::cout << float(TCOUNT) / float (STAT) << std::endl;
	
    return 0;
	
}
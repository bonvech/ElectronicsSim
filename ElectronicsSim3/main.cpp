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

/**
 * @brief Класс для моделирования электроники
 *
 * Данный класс используется для загрузки входных параметров
 * и расчета конечных данных.
 *
 */
class ModelElectronics {
private:
    /* Константы для расчетов. */
    const int N_CHAN = 2653; /// Число каналов
    const int PULSE_LENGTH = 10000;    
    const int BIN_2_GEN = 500;
	const int TicksPerNs = 10;    /// сколько тиков в наносекунду гененрировать (частота моделирования в ГГц)
	const int BinLength = 10;     /// Сколько наносекунд в одном бине оцифровки
    //const float PH_Flux = 1.222;  /// Ожидаемый поток фоновых фотонов в 1 нс ( 0.08 вычислен из потока для СФЕРЫ-2 с учётом площади диафрагмы, угла обзора и квантовой эффективности; по Галкину 1.222)
	const float TargetAmp = 7;   /// Целевая амплитуда однофотоэлектронного пика 7 ADU (при температуре -15) 
	const float TargetUover = 5; /// Целевое перенапряжение 5В при температуре -15 градусов
    std::vector<float> curbase; /// Относительные токи
    std::vector<float> pulse; /// Импульсные характеристики тока    
	std::vector<float> params; /// Строка парамтеров состояния детектора
    std::vector<float> Toff; /// Относительные сдвиги каналов            
    int N_PHEL{}; /// Число зарегистрированных фотонов
    std::vector<int> PMTid; /// Номера ФЭУ, куда упали фотоэлектроны
    std::vector<float> T; /// Времена прихода фотоэлектронов
	std::vector<int> PhType; /// Времена прихода фотоэлектронов
    float Tmin{}; /// Минимальное время прихода    
    std::vector<std::vector<int>> data_out; /// Выводной массив
    std::vector<std::vector<float>> data; /// Данные
	float U; /// напражение на SiPM
	float Temp; /// температура SiPM
	float cross_exp;
	float CurAmp;
	float PH_Flux;
public:
    ModelElectronics();

    void GetMoshits();

    void GetSimple(const std::string &, std::vector<float> &);

//    void GetC();

	void SetUpDetector();

    void GenerateEvent();

    void AddBackground();

    void SimulateDig();

    void PrintDataOut();

    std::vector<float> &getCurbaseRef() { return curbase; }

    std::vector<float> &getPulseRef() { return pulse; }
	
	std::vector<float> &getToffRef() { return Toff; }
	
	std::vector<float> &getParamsRef() { return params; }
};

ModelElectronics::ModelElectronics() :        
        curbase(N_CHAN, 0),
        pulse(PULSE_LENGTH, 0),         
        Toff(N_CHAN, 0),
		params(4,0),
        data_out(N_CHAN, std::vector<int>(BIN_2_GEN, 0)),
        data(N_CHAN, std::vector<float>(BIN_2_GEN * TicksPerNs * BinLength + 2 * PULSE_LENGTH + 2, 0)) {}
/**
 * @brief Функция для считывания файла moshits
 */
void ModelElectronics::GetMoshits() {
    std::ifstream moshits("mosaic_hits");

    if (!moshits.is_open()) {
        std::cerr << "Failed to open the moshits file!" << std::endl;
    }
    N_PHEL = int(std::count(std::istreambuf_iterator<char>(moshits),
                            std::istreambuf_iterator<char>(), '\n')) - 1;

    moshits.clear();
    moshits.seekg(0, std::ios::beg);

    std::string line;
    std::getline(moshits, line);
    double m_tmp, s_tmp, t_tmp, b_tmp, tmp;
    while (std::getline(moshits, line)) {
        std::istringstream ss(line);
        ss >> m_tmp;
		ss >> s_tmp;
        for (int i{0}; i < 3; i++) {
            ss >> tmp;
        }
        ss >> t_tmp;
		for (int i{0}; i < 3; i++) {
            ss >> tmp;
        }
		ss >> b_tmp;
		PhType.push_back(int(b_tmp));
		PMTid.push_back(int(m_tmp)*7+int(s_tmp));
		T.push_back(float(t_tmp));		
    }
    auto min_it{std::min_element(T.begin(), T.end())};
    Tmin = *min_it;
    moshits.close();
}

/**
 * @brief Обобщенная функция для считывания файлов с одной строкой
 * @param fileName Имя файла
 * @param vec ссылка на массив, в который надо считать файл
 */
void ModelElectronics::GetSimple(const std::string &fileName, std::vector<float> &vec) {
    std::ifstream input(fileName);
    if (!input.is_open()) {
        std::cerr << "Failed to open the " << fileName << " file!" << std::endl;
    }
    std::string line;
    std::getline(input, line);
    std::istringstream ss(line);
    auto it{vec.begin()};
    float tmp;
    while (ss >> tmp && it != vec.end()) {
        *it = tmp;
        ++it;
    }
    input.close();
}

/**
 * @brief Функция для считывания файла калибровки
 
void ModelElectronics::GetC() {
    std::ifstream cal("14484.cal");
    if (!cal.is_open()) {
        std::cerr << "Failed to open the calibration file!" << std::endl;
    }
    std::string line;
    float tmp;
    auto it{Cal.begin()};
    while (std::getline(cal, line) && it != Cal.end()) {
        std::istringstream ss(line);
        for (int i{0}; i < 13; i++) {
            ss >> tmp;
        }
        ss >> *it;
        ++it;
    }
    cal.close();
}
*/


/**
 * @brief Настройка матрицы, вероятность кросстока, коэффициент усиления
 */ 
void ModelElectronics::SetUpDetector() {
	float Ubr;
	float crossprob;
	// файл params по порядку содержит
	// 1. Напряжение на SiPM в Вольтах
	// 2. Температуру мозаики в градусах
	// 3. Ширину распределения коэффициента усилиения в условных единцах (надо уточнять, но оно около, хотя может быть будет куда как лушче 0.8)
	// 4. Параметры фона в единицах фотон на пиксель за 1 нс
	
	U = params[0];   // чтение напряжения из параметров	
	Temp = params[1]; // чтение температуры из параметров
	Ubr = 23.9 + Temp * 21.5e-3;  // перенапряжение пробоя (не совсем паспорт SiPM, скорее диплом Аминевой А.А.)
	crossprob = (Ubr-2.0) * 0.1; // вероятность кросстоков (см. диплом Аминевой А.А.)	
	crossprob = crossprob * 0.7; // кросс-токи для SiPM без светосборника будут ощутимо выше, чем они же для ситуации со светосборником (в силу физики процесса)
	if (crossprob < 0.04) {  
		crossprob = 0.04;  // так в паспорте
	}
	// вероятность кросс-токов - это вероятность получить не нулевое количество дополнительных сработавших 
	// ячеек при срабатывании текущей, т.е. по Пуассону P(0)=1-p=exp(-lambda)
	// lambda = - ln (1-p)
	
	cross_exp = log (1 - crossprob);	
	CurAmp = TargetAmp * ( U - Ubr) / TargetUover;
	
	PH_Flux=params[3];  // поток фотонов в пиксель
}

/**
 * @brief Генерация события
 */
 
void ModelElectronics::GenerateEvent() {
    int T_ph;
	int N;
    float amp_ph;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<> dist(cross_exp);
	std::normal_distribution<> subamp(0,params[2]);
	
    for (int phid{0}; phid < N_PHEL; phid++) {
		if (PhType[phid] > 1) {
			continue;
		}
		N = 1;
		amp_ph = 0;
		while (N > 0) {
			N += (dist(gen) - 1);
			amp_ph += 1;			
		}
        amp_ph += subamp(gen);		
        T_ph = int(TicksPerNs * (T[phid] - Tmin) + floor(0.45 * BIN_2_GEN * TicksPerNs * BinLength + PULSE_LENGTH));
        for (int t{0}; t < PULSE_LENGTH; t++) {
            data[PMTid[phid]][t + T_ph] += CurAmp * amp_ph * pulse[t];
        }
    }
}

/**
 * @brief Добавление фона
 */
void ModelElectronics::AddBackground() {
    int BG_LENGTH{ (BIN_2_GEN + 1) * TicksPerNs * BinLength + PULSE_LENGTH};
    float N_AVG;
    float N_PHEL_exp;
    float amp_ph;
    int T_ph;
	int N;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<> dist(cross_exp);
	std::normal_distribution<> subamp(0,params[2]);
    std::uniform_int_distribution<> dist_bg(0, BG_LENGTH);
	
    for (int j{0}; j < N_CHAN; j++) {		
        N_AVG = PH_Flux * curbase[j];		
        N_PHEL_exp = N_AVG * float(BG_LENGTH) / TicksPerNs;		
        std::poisson_distribution bgphdist(N_PHEL_exp);
        N_PHEL = bgphdist(gen);		
        for (int n{0}; n < N_PHEL; n++) {
            
			N = 1;
			amp_ph = 0;
			while (N > 0) {
				N += (dist(gen) - 1);
				amp_ph += 1;			
			}
			amp_ph += subamp(gen);	
			
            T_ph = dist_bg(gen);
            for (int t{0}; t < PULSE_LENGTH; t++) {
                data[j][t + T_ph] += CurAmp * amp_ph * pulse[t];
            }
        }
    }
}

/**
 * @brief Имитация оцифровки
 */
void ModelElectronics::SimulateDig() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist_shift(0, TicksPerNs * BinLength);
    int t_shift{dist_shift(gen)};  
    for (int j{0}; j < N_CHAN; j++) {        
        t_shift = dist_shift(gen);     
        for (int i{0}; i < BIN_2_GEN; i++) {
            data_out[j][i] = int(data[j][PULSE_LENGTH + t_shift + int (i + Toff[j]) * TicksPerNs * BinLength ]);
        }
    }

}

/**
 * @brief Метод для печати выводного файла
 */
void ModelElectronics::PrintDataOut() {
    std::ofstream outFile("data_out");
    if (!outFile.is_open()) {
        std::cerr << "Open file error." << std::endl;
    }
    for (const auto& innerVec : data_out){
        for (const auto& item: innerVec){
            outFile << item << ' ';
        }
        outFile << std::endl;
    }
    outFile.close();
}

/**
 * @brief Класс для многопоточности
 */
class ThreadManager {
private:
    ModelElectronics &obj;
public:
    explicit ThreadManager(ModelElectronics &objRef) : obj(objRef) {}

    void inputAll();
};

/**
 * @brief Метод для заполнения массивов в многопоточном режиме
 */
void ThreadManager::inputAll() {
    std::vector<std::thread> threads;
//    threads.emplace_back(&ModelElectronics::GetC, &obj);
    threads.emplace_back(&ModelElectronics::GetMoshits, &obj);
    threads.emplace_back(&ModelElectronics::GetSimple, &obj, "CurRels.dat", std::ref(obj.getCurbaseRef()));
    threads.emplace_back(&ModelElectronics::GetSimple, &obj, "Impulse10GHz.dat", std::ref(obj.getPulseRef()));
	threads.emplace_back(&ModelElectronics::GetSimple, &obj, "Toff.dat", std::ref(obj.getToffRef()));
	threads.emplace_back(&ModelElectronics::GetSimple, &obj, "Params.dat", std::ref(obj.getParamsRef()));
    for (auto &t: threads) {
        t.join();
    }
}


int main() {
    ModelElectronics model;	
    ThreadManager manager(model);;
    manager.inputAll();
	model.SetUpDetector();
    model.GenerateEvent();
    model.AddBackground();	
    model.SimulateDig();
    model.PrintDataOut();
    return 0;
}
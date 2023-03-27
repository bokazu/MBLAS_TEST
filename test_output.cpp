#include "test_output.hpp"

using namespace std;

void test_mm_ddot(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 初期状態行列の用意
    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.001);
    uniform_real_distribution<> rand2(0, 0.001);

    double **M1 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M1[i] = new double[mat_dim];

    double **M2 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M2[i] = new double[mat_dim];

    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M1[i][j] = rand1(mt);
            M2[i][j] = rand2(mt);
        }
    }

    // 内積計算の実施と計算時間の測定を行う
    chrono::system_clock::time_point mm_ddot_start, MP_mm_ddot_start, MP_schedule_mm_ddot_start, mm_ddot_end, MP_mm_ddot_end, MP_schedule_mm_ddot_end;
    cout << "\n\n Test : mm_ddot()\n";
    cout << "====================================================\n";
    // 各種内積計算用の計算結果を格納する
    vector<double> res_mm_ddot;
    vector<double> res_MP_ddot;
    vector<double> res_MP_schedule_ddot;
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mm_ddot;
    vector<double> time_MP_ddot;
    vector<double> time_MP_schedule_ddot;

    // trial_num回だけテストを実施する
    for (int i = 0; i < trial_num; i++)
    {
        // About mm_ddot
        mm_ddot_start = chrono::system_clock::now();
        double mat_ddot = mm_ddot(mat_dim, M1, M2);
        mm_ddot_end = chrono::system_clock::now();
        // About MP_mm_ddot
        MP_mm_ddot_start = chrono::system_clock::now();
        double MP_mat_ddot = MP_mm_ddot(mat_dim, M1, M2);
        MP_mm_ddot_end = chrono::system_clock::now();
        // About MP_schedule_mm_ddot
        MP_schedule_mm_ddot_start = chrono::system_clock::now();
        double MP_schedule_mat_ddot = MP_schedule_mm_ddot(mat_dim, M1, M2);
        MP_schedule_mm_ddot_end = chrono::system_clock::now();

        // 時間計測結果をミリ秒に変換
        auto mm_ddot_time_msec = chrono::duration_cast<chrono::milliseconds>(mm_ddot_end - mm_ddot_start).count();
        auto MP_mm_ddot_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_mm_ddot_end - MP_mm_ddot_start).count();
        auto MP_schedule_mm_ddot_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_schedule_mm_ddot_end - MP_schedule_mm_ddot_start).count();

        /*計算結果を各種配列に保存する*/
        res_mm_ddot.push_back(mat_ddot);
        res_MP_ddot.push_back(MP_mat_ddot);
        res_MP_schedule_ddot.push_back(MP_schedule_mat_ddot);
        /*時間測定結果を各種配列に保存する*/
        time_mm_ddot.push_back(mm_ddot_time_msec);
        time_MP_ddot.push_back(MP_mm_ddot_time_msec);
        time_MP_schedule_ddot.push_back(MP_schedule_mm_ddot_time_msec);
    }

    // 計算結果をファイル出力する
    std::ofstream ofs_calc(PATH_result_of_calc);
    // fileをopenできるか確認する
    if (!ofs_calc)
    {
        printf("@test1_mm_ddot error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << res_mm_ddot[i] << "," << res_MP_ddot[i] << "," << res_MP_schedule_ddot[i] << std::endl;
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << res_mm_ddot[i] << " " << res_MP_ddot[i] << " " << res_MP_schedule_ddot[i] << std::endl;
            }
        }
    }
    ofs_calc.close();

    // 時間測定結果をファイル出力する
    std::ofstream ofs_time(PATH_result_of_time);
    // fileをopenできるか確認する
    if (!ofs_time)
    {
        printf("@test1_mm_ddot error:: \"%s\" could not open.", PATH_result_of_time.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            ofs_time << "Normal"
                     << ","
                     << "OpenMP"
                     << "OpenMP with schedule" << endl;
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_ddot[i] << "," << time_MP_ddot[i] << "," << time_MP_schedule_ddot[i] << std::endl;
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_ddot[i] << " " << time_MP_ddot[i] << " " << time_MP_schedule_ddot[i] << std::endl;
            }
        }
    }
    ofs_time.close();

    /*---------------------- -メモリの開放-----------------------*/
    for (int i = 0; i < mat_dim; i++)
        delete[] M1[i];
    delete[] M1;

    for (int i = 0; i < mat_dim; i++)
        delete[] M2[i];
    delete[] M2;
};

void test_mm_dcopy(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 初期状態行列の用意
    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(-0.5, 0.5);

    double **M1 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M1[i] = new double[mat_dim];

    double **M2 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M2[i] = new double[mat_dim];

    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M1[i][j] = rand1(mt);
            M2[i][j] = 0.;
        }
    }

    // 内積計算の実施と計算時間の測定を行う
    chrono::system_clock::time_point mm_dcopy_start, MP_mm_dcopy_start, MP_schedule_mm_dcopy_start, mm_dcopy_end, MP_mm_dcopy_end, MP_schedule_mm_dcopy_end;
    cout << "\n\n Test : mm_dcopy()\n";
    cout << "====================================================\n";
    // 各種内積計算用の計算結果を格納する
    bool is_ok_mm_dcopy = true;
    bool is_ok_MP_mm_dcopy = true;
    bool is_ok_MP_schedule_mm_dcopy = true;
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mm_dcopy;
    vector<double> time_MP_dcopy;
    vector<double> time_MP_schedule_dcopy;

    // trial_num回だけテストを実施する
    for (int i = 0; i < trial_num; i++)
    {
        // About mm_dcopy
        mm_dcopy_start = chrono::system_clock::now();
        mm_dcopy(mat_dim, M1, M2);
        mm_dcopy_end = chrono::system_clock::now();

        // 正しくコピーを行うことができているかの確認
        for (int i = 0; i < mat_dim; i++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                if (M1[i][j] != M2[i][j])
                    is_ok_mm_dcopy = false;
            }
        }
// 行列M2の初期化
#pragma omp parallel for
        for (int i = 0; i < mat_dim; i++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                M2[i][j] = 0.0;
            }
        }
        // 正しくコピーを行うことができているかのチェック

        //  About MP_mm_dcopy
        MP_mm_dcopy_start = chrono::system_clock::now();
        MP_mm_dcopy(mat_dim, M1, M2);
        MP_mm_dcopy_end = chrono::system_clock::now();

        // 正しくコピーを行うことができているかの確認
        for (int i = 0; i < mat_dim; i++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                if (M1[i][j] != M2[i][j])
                    is_ok_MP_mm_dcopy = false;
            }
        }
        // 行列M2の初期化
#pragma omp parallel for
        for (int i = 0; i < mat_dim; i++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                M2[i][j] = 0.0;
            }
        }
        // About MP_schedule_mm_dcopy
        MP_schedule_mm_dcopy_start = chrono::system_clock::now();
        MP_schedule_mm_dcopy(mat_dim, M1, M2);
        MP_schedule_mm_dcopy_end = chrono::system_clock::now();

        // 正しくコピーを行うことができているかの確認
        for (int i = 0; i < mat_dim; i++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                if (M1[i][j] != M2[i][j])
                    is_ok_MP_schedule_mm_dcopy = false;
            }
        }

        // 時間計測結果をミリ秒に変換
        auto mm_dcopy_time_msec = chrono::duration_cast<chrono::milliseconds>(mm_dcopy_end - mm_dcopy_start).count();
        auto MP_mm_dcopy_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_mm_dcopy_end - MP_mm_dcopy_start).count();
        auto MP_schedule_mm_dcopy_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_schedule_mm_dcopy_end - MP_schedule_mm_dcopy_start).count();

        /*計算結果を各種配列に保存する*/

        /*時間測定結果を各種配列に保存する*/
        time_mm_dcopy.push_back(mm_dcopy_time_msec);
        time_MP_dcopy.push_back(MP_mm_dcopy_time_msec);
        time_MP_schedule_dcopy.push_back(MP_schedule_mm_dcopy_time_msec);
    }

    // 計算結果をファイル出力する
    std::ofstream ofs_calc(PATH_result_of_calc);
    // fileをopenできるか確認する
    if (!ofs_calc)
    {
        printf("@test1_mm_dcopy error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                if (is_ok_mm_dcopy)
                    ofs_calc << "success"
                             << ",";
                else
                    ofs_calc << "failed"
                             << ",";

                if (is_ok_MP_mm_dcopy)
                    ofs_calc << "success"
                             << ",";
                else
                    ofs_calc << "failed"
                             << ",";

                if (is_ok_MP_schedule_mm_dcopy)
                    ofs_calc << "success\n";
                else
                    ofs_calc << "failed\n";
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                if (is_ok_mm_dcopy)
                    ofs_calc << "success"
                             << " ";
                else
                    ofs_calc << "failed"
                             << " ";

                if (is_ok_MP_mm_dcopy)
                    ofs_calc << "success"
                             << " ";
                else
                    ofs_calc << "failed"
                             << " ";

                if (is_ok_MP_schedule_mm_dcopy)
                    ofs_calc << "success\n";
                else
                    ofs_calc << "failed\n";
            }
        }
    }
    ofs_calc.close();

    // 時間測定結果をファイル出力する
    std::ofstream ofs_time(PATH_result_of_time);
    // fileをopenできるか確認する
    if (!ofs_time)
    {
        printf("@test1_mm_dcopy error:: \"%s\" could not open.", PATH_result_of_time.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            ofs_time << "Normal"
                     << ","
                     << "OpenMP"
                     << "OpenMP with schedule" << endl;
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_dcopy[i] << "," << time_MP_dcopy[i] << "," << time_MP_schedule_dcopy[i] << std::endl;
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_dcopy[i] << " " << time_MP_dcopy[i] << " " << time_MP_schedule_dcopy[i] << std::endl;
            }
        }
    }
    ofs_time.close();

    /*---------------------- -メモリの開放-----------------------*/
    for (int i = 0; i < mat_dim; i++)
        delete[] M1[i];
    delete[] M1;

    for (int i = 0; i < mat_dim; i++)
        delete[] M2[i];
    delete[] M2;
};

void test_mm_dnrm2(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 初期状態行列の用意
    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.001);
    uniform_real_distribution<> rand2(0, 0.001);

    double **M1 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M1[i] = new double[mat_dim];

#pragma omp parallel for
    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M1[i][j] = rand1(mt);
        }
    }

    // 内積計算の実施と計算時間の測定を行う
    chrono::system_clock::time_point mm_dnrm2_start, MP_mm_dnrm2_start, MP_schedule_mm_dnrm2_start, mm_dnrm2_end, MP_mm_dnrm2_end, MP_schedule_mm_dnrm2_end;
    cout << "\n\n Test : mm_dnrm2()\n";
    cout << "====================================================\n";
    // 各種内積計算用の計算結果を格納する
    vector<double> res_mm_dnrm2;
    vector<double> res_MP_dnrm2;
    vector<double> res_MP_schedule_dnrm2;
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mm_dnrm2;
    vector<double> time_MP_dnrm2;
    vector<double> time_MP_schedule_dnrm2;

    // trial_num回だけテストを実施する
    for (int i = 0; i < trial_num; i++)
    {
        // About mm_dnrm2
        mm_dnrm2_start = chrono::system_clock::now();
        double mat_dnrm2 = mm_dnrm2(mat_dim, M1);
        mm_dnrm2_end = chrono::system_clock::now();
        // About MP_mm_dnrm2
        MP_mm_dnrm2_start = chrono::system_clock::now();
        double MP_mat_dnrm2 = MP_mm_dnrm2(mat_dim, M1);
        MP_mm_dnrm2_end = chrono::system_clock::now();
        // About MP_schedule_mm_dnrm2
        MP_schedule_mm_dnrm2_start = chrono::system_clock::now();
        double MP_schedule_mat_dnrm2 = MP_schedule_mm_dnrm2(mat_dim, M1);
        MP_schedule_mm_dnrm2_end = chrono::system_clock::now();

        // 時間計測結果をミリ秒に変換
        auto mm_dnrm2_time_msec = chrono::duration_cast<chrono::milliseconds>(mm_dnrm2_end - mm_dnrm2_start).count();
        auto MP_mm_dnrm2_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_mm_dnrm2_end - MP_mm_dnrm2_start).count();
        auto MP_schedule_mm_dnrm2_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_schedule_mm_dnrm2_end - MP_schedule_mm_dnrm2_start).count();

        /*計算結果を各種配列に保存する*/
        res_mm_dnrm2.push_back(mat_dnrm2);
        res_MP_dnrm2.push_back(MP_mat_dnrm2);
        res_MP_schedule_dnrm2.push_back(MP_schedule_mat_dnrm2);
        /*時間測定結果を各種配列に保存する*/
        time_mm_dnrm2.push_back(mm_dnrm2_time_msec);
        time_MP_dnrm2.push_back(MP_mm_dnrm2_time_msec);
        time_MP_schedule_dnrm2.push_back(MP_schedule_mm_dnrm2_time_msec);
    }

    // 計算結果をファイル出力する
    std::ofstream ofs_calc(PATH_result_of_calc);
    // fileをopenできるか確認する
    if (!ofs_calc)
    {
        printf("@test1_mm_dnrm2 error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << res_mm_dnrm2[i] << "," << res_MP_dnrm2[i] << "," << res_MP_schedule_dnrm2[i] << std::endl;
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << res_mm_dnrm2[i] << " " << res_MP_dnrm2[i] << " " << res_MP_schedule_dnrm2[i] << std::endl;
            }
        }
    }
    ofs_calc.close();

    // 時間測定結果をファイル出力する
    std::ofstream ofs_time(PATH_result_of_time);
    // fileをopenできるか確認する
    if (!ofs_time)
    {
        printf("@test1_mm_dnrm2 error:: \"%s\" could not open.", PATH_result_of_time.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            ofs_time << "Normal"
                     << ","
                     << "OpenMP"
                     << "OpenMP with schedule" << endl;
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_dnrm2[i] << "," << time_MP_dnrm2[i] << "," << time_MP_schedule_dnrm2[i] << std::endl;
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_dnrm2[i] << " " << time_MP_dnrm2[i] << " " << time_MP_schedule_dnrm2[i] << std::endl;
            }
        }
    }
    ofs_time.close();

    /*---------------------- -メモリの開放-----------------------*/
    for (int i = 0; i < mat_dim; i++)
        delete[] M1[i];
    delete[] M1;
};

void test_mm_dscal(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{

    // 初期状態行列の用意
    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.001);

    double **M1 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M1[i] = new double[mat_dim];

    double **M2 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M2[i] = new double[mat_dim];

    double **M3 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M3[i] = new double[mat_dim];

    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M1[i][j] = rand1(mt);
            M2[i][j] = M1[i][j];
            M3[i][j] = M1[i][j];
        }
    }

    // 内積計算の実施と計算時間の測定を行う
    chrono::system_clock::time_point mm_dscal_start, MP_mm_dscal_start, MP_schedule_mm_dscal_start, mm_dscal_end, MP_mm_dscal_end, MP_schedule_mm_dscal_end;
    cout << "\n\n Test : mm_dscal()\n";
    cout << "====================================================\n";
    // 各種内積計算用の計算結果を格納する
    bool is_ok_mm_dscal = true;
    bool is_ok_MP_mm_dscal = true;
    bool is_ok_MP_schedule_mm_dscal = true;
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mm_dscal;
    vector<double> time_MP_dscal;
    vector<double> time_MP_schedule_dscal;

    // trial_num回だけテストを実施する
    for (int i = 0; i < trial_num; i++)
    {
        //[3/10]計算結果の確認方法を修正する
        // About mm_dscal
        mm_dscal_start = chrono::system_clock::now();
        mm_dscal(mat_dim, 0.5, M1);
        mm_dscal_end = chrono::system_clock::now();

        // //  About MP_mm_dscal
        MP_mm_dscal_start = chrono::system_clock::now();
        MP_mm_dscal(mat_dim, 0.5, M2);
        MP_mm_dscal_end = chrono::system_clock::now();

        // // 正しく並列計算を行うことができているかの確認
        cout << "Now checking result of calculation. Please wait ...\n";
        for (int i = 0; i < mat_dim; i++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                if (M1[i][j] != M2[i][j])
                    is_ok_MP_mm_dscal = false;
            }
        }

        if (is_ok_mm_dscal)
            cout << "Success!\n";
        else
            cout << "Failed.\n";

        // // About MP_schedule_mm_dscal
        MP_schedule_mm_dscal_start = chrono::system_clock::now();
        MP_schedule_mm_dscal(mat_dim, 0.5, M3);
        MP_schedule_mm_dscal_end = chrono::system_clock::now();

        // // 正しく並列計算を行うことができているかの確認
        cout << "Now checking result of calculation. Please wait ...\n";
        for (int i = 0; i < mat_dim; i++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                if (M1[i][j] != M3[i][j])
                    is_ok_MP_schedule_mm_dscal = false;
            }
        }
        if (is_ok_MP_mm_dscal)
            cout << "Success!\n";
        else
            cout << "Failed.\n";

        // 時間計測結果をミリ秒に変換
        auto mm_dscal_time_msec = chrono::duration_cast<chrono::milliseconds>(mm_dscal_end - mm_dscal_start).count();
        auto MP_mm_dscal_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_mm_dscal_end - MP_mm_dscal_start).count();
        auto MP_schedule_mm_dscal_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_schedule_mm_dscal_end - MP_schedule_mm_dscal_start).count();

        /*計算結果を各種配列に保存する*/

        /*時間測定結果を各種配列に保存する*/
        time_mm_dscal.push_back(mm_dscal_time_msec);
        time_MP_dscal.push_back(MP_mm_dscal_time_msec);
        time_MP_schedule_dscal.push_back(MP_schedule_mm_dscal_time_msec);
    }

    // 計算結果をファイル出力する
    std::ofstream ofs_calc(PATH_result_of_calc);
    // fileをopenできるか確認する
    if (!ofs_calc)
    {
        printf("@test1_mm_dscal error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                if (is_ok_mm_dscal)
                    ofs_calc << "success"
                             << ",";
                else
                    ofs_calc << "failed"
                             << ",";

                if (is_ok_MP_mm_dscal)
                    ofs_calc << "success"
                             << ",";
                else
                    ofs_calc << "failed"
                             << ",";

                if (is_ok_MP_schedule_mm_dscal)
                    ofs_calc << "success\n";
                else
                    ofs_calc << "failed\n";
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                if (is_ok_mm_dscal)
                    ofs_calc << "success"
                             << " ";
                else
                    ofs_calc << "failed"
                             << " ";

                if (is_ok_MP_mm_dscal)
                    ofs_calc << "success"
                             << " ";
                else
                    ofs_calc << "failed"
                             << " ";

                if (is_ok_MP_schedule_mm_dscal)
                    ofs_calc << "success\n";
                else
                    ofs_calc << "failed\n";
            }
        }
    }
    ofs_calc.close();

    // 時間測定結果をファイル出力する
    std::ofstream ofs_time(PATH_result_of_time);
    // fileをopenできるか確認する
    if (!ofs_time)
    {
        printf("@test1_mm_dscal error:: \"%s\" could not open.", PATH_result_of_time.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            ofs_time << "Normal"
                     << ","
                     << "OpenMP"
                     << "OpenMP with schedule" << endl;
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_dscal[i] << "," << time_MP_dscal[i] << "," << time_MP_schedule_dscal[i] << std::endl;
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_dscal[i] << " " << time_MP_dscal[i] << " " << time_MP_schedule_dscal[i] << std::endl;
            }
        }
    }
    ofs_time.close();

    /*---------------------- -メモリの開放-----------------------*/
    for (int i = 0; i < mat_dim; i++)
        delete[] M1[i];
    delete[] M1;

    for (int i = 0; i < mat_dim; i++)
        delete[] M2[i];
    delete[] M2;

    for (int i = 0; i < mat_dim; i++)
        delete[] M3[i];
    delete[] M3;
};

void test_mm_daxpy(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 初期状態行列の用意
    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.001);

    double **M1 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M1[i] = new double[mat_dim];

    double **M2 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M2[i] = new double[mat_dim];

    double **M3 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M3[i] = new double[mat_dim];

    double **M4 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M4[i] = new double[mat_dim];

    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M1[i][j] = rand1(mt);
            M2[i][j] = 0.0;
            M3[i][j] = 0.0;
            M4[i][j] = 0.0;
        }
    }

    // daxpyの実施と計算時間の測定を行う
    chrono::system_clock::time_point mm_daxpy_start, MP_mm_daxpy_start, MP_schedule_mm_daxpy_start, mm_daxpy_end, MP_mm_daxpy_end, MP_schedule_mm_daxpy_end;
    cout << "\n\n Test : mm_daxpy()\n";
    cout << "====================================================\n";
    // daxpyの計算結果について適当に行列の（(1,3)成分を選び結果の確認に使用する
    vector<double> M2_13;
    vector<double> M3_13;
    vector<double> M4_13;
    // daxpyの計算結果について適当に行列の（(10,10)成分を選び結果の確認に使用する
    vector<double> M2_tt;
    vector<double> M3_tt;
    vector<double> M4_tt;
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mm_daxpy;
    vector<double> time_MP_daxpy;
    vector<double> time_MP_schedule_daxpy;

    // // trial_num回だけテストを実施する
    for (int i = 0; i < trial_num; i++)
    {
        // About mm_daxpy
        mm_daxpy_start = chrono::system_clock::now();
        mm_daxpy(mat_dim, 0.5, M1, M2);
        mm_daxpy_end = chrono::system_clock::now();
        M2_13.push_back(M2[1][3]);
        M2_tt.push_back(M2[10][10]);

        //  About MP_mm_daxpy
        MP_mm_daxpy_start = chrono::system_clock::now();
        MP_mm_daxpy(mat_dim, 0.5, M1, M3);
        MP_mm_daxpy_end = chrono::system_clock::now();
        M3_13.push_back(M3[1][3]);
        M3_tt.push_back(M3[10][10]);

        // About MP_schedule_mm_daxpy
        MP_schedule_mm_daxpy_start = chrono::system_clock::now();
        MP_schedule_mm_daxpy(mat_dim, 0.5, M1, M4);
        MP_schedule_mm_daxpy_end = chrono::system_clock::now();
        M4_13.push_back(M4[1][3]);
        M4_tt.push_back(M4[10][10]);

        // 時間計測結果をミリ秒に変換
        auto mm_daxpy_time_msec = chrono::duration_cast<chrono::milliseconds>(mm_daxpy_end - mm_daxpy_start).count();
        auto MP_mm_daxpy_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_mm_daxpy_end - MP_mm_daxpy_start).count();
        auto MP_schedule_mm_daxpy_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_schedule_mm_daxpy_end - MP_schedule_mm_daxpy_start).count();

        /*計算結果を各種配列に保存する*/

        /*時間測定結果を各種配列に保存する*/
        time_mm_daxpy.push_back(mm_daxpy_time_msec);
        time_MP_daxpy.push_back(MP_mm_daxpy_time_msec);
        time_MP_schedule_daxpy.push_back(MP_schedule_mm_daxpy_time_msec);
    }

    // 計算結果をファイル出力する
    std::ofstream ofs_calc(PATH_result_of_calc);
    // fileをopenできるか確認する
    if (!ofs_calc)
    {
        printf("@test1_mm_dnrm2 error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            ofs_calc << "------------------------------------------------------------------------------\n";
            ofs_calc << "@M[1][3]\n";
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << M2_13[i] << "," << M3_13[i] << "," << M4_13[i] << std::endl;
            }
            ofs_calc << "------------------------------------------------------------------------------\n";

            ofs_calc << "------------------------------------------------------------------------------\n";
            ofs_calc << "@M[10][10]\n";
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << M2_tt[i] << "," << M3_tt[i] << "," << M4_tt[i] << std::endl;
            }
            ofs_calc << "------------------------------------------------------------------------------\n";
        }
        else // txt出力
        {
            ofs_calc << "------------------------------------------------------------------------------\n";
            ofs_calc << "@M[1][3]\n";
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << M2_13[i] << " " << M3_13[i] << " " << M4_13[i] << std::endl;
            }
            ofs_calc << "------------------------------------------------------------------------------\n";

            ofs_calc << "------------------------------------------------------------------------------\n";
            ofs_calc << "@M[10][10]\n";
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << M2_tt[i] << " " << M3_tt[i] << " " << M4_tt[i] << std::endl;
            }
            ofs_calc << "------------------------------------------------------------------------------\n";
        }
    }
    ofs_calc.close();

    // 時間測定結果をファイル出力する
    std::ofstream ofs_time(PATH_result_of_time);
    // fileをopenできるか確認する
    if (!ofs_time)
    {
        printf("@test1_mm_dnrm2 error:: \"%s\" could not open.", PATH_result_of_time.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            ofs_time << "Normal"
                     << ","
                     << "OpenMP"
                     << "OpenMP with schedule" << endl;
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_daxpy[i] << "," << time_MP_daxpy[i] << "," << time_MP_schedule_daxpy[i] << std::endl;
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_daxpy[i] << " " << time_MP_daxpy[i] << " " << time_MP_schedule_daxpy[i] << std::endl;
            }
        }
    }
    ofs_time.close();

    /*---------------------- -メモリの開放-----------------------*/
    for (int i = 0; i < mat_dim; i++)
        delete[] M1[i];
    delete[] M1;

    for (int i = 0; i < mat_dim; i++)
        delete[] M2[i];
    delete[] M2;

    for (int i = 0; i < mat_dim; i++)
        delete[] M3[i];
    delete[] M3;

    for (int i = 0; i < mat_dim; i++)
        delete[] M4[i];
    delete[] M4;
};

void test_mm_sdz(int trial_num, int mat_dim, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 初期状態行列の用意
    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.001);

    // 規格化する状態行列
    double **M1 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M1[i] = new double[mat_dim];

    // 規格化する状態行列1
    double **M2 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M2[i] = new double[mat_dim];

    // 規格化する状態行列2
    double **M3 = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
        M3[i] = new double[mat_dim];

    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = 0; j < mat_dim; j++)
        {
            M1[i][j] = rand1(mt);
            M2[i][j] = M1[i][j];
            M3[i][j] = M1[i][j];
        }
    }

    // daxpyの実施と計算時間の測定を行う
    chrono::system_clock::time_point mm_sdz_start, MP_mm_sdz_start, MP_schedule_mm_sdz_start, mm_sdz_end, MP_mm_sdz_end, MP_schedule_mm_sdz_end;
    cout << "\n\n Test : mm_sdz()\n";
    cout << "====================================================\n";
    // sdzの計算結果について適当に状態行列の（(1,3)成分を選び結果の確認に使用する
    vector<double> M2_13;
    vector<double> M3_13;
    vector<double> M4_13;
    // sdzの計算結果について適当に状態行列の（(10,10)成分を選び結果の確認に使用する
    vector<double> M2_tt;
    vector<double> M3_tt;
    vector<double> M4_tt;
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mm_sdz;
    vector<double> time_MP_sdz;
    vector<double> time_MP_schedule_sdz;

    // // trial_num回だけテストを実施する
    for (int i = 0; i < trial_num; i++)
    {
        // About mm_sdz
        mm_sdz_start = chrono::system_clock::now();
        mm_sdz(mat_dim, M1);
        mm_sdz_end = chrono::system_clock::now();
        M2_13.push_back(M1[1][3]);
        M2_tt.push_back(M1[10][10]);

        // M1をもとの行列に初期化する
        MP_mm_dcopy(mat_dim, M2, M1);

        //  About MP_mm_sdz
        MP_mm_sdz_start = chrono::system_clock::now();
        MP_mm_sdz(mat_dim, M2);
        MP_mm_sdz_end = chrono::system_clock::now();
        M3_13.push_back(M2[1][3]);
        M3_tt.push_back(M2[10][10]);

        // M2をもとの行列に初期化する
        MP_mm_dcopy(mat_dim, M1, M2);

        // About MP_schedule_mm_sdz
        MP_schedule_mm_sdz_start = chrono::system_clock::now();
        MP_schedule_mm_sdz(mat_dim, M3);
        MP_schedule_mm_sdz_end = chrono::system_clock::now();
        M4_13.push_back(M3[1][3]);
        M4_tt.push_back(M3[10][10]);

        // M2をもとの行列に初期化する
        MP_mm_dcopy(mat_dim, M1, M3);

        // 時間計測結果をミリ秒に変換
        auto mm_sdz_time_msec = chrono::duration_cast<chrono::milliseconds>(mm_sdz_end - mm_sdz_start).count();
        auto MP_mm_sdz_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_mm_sdz_end - MP_mm_sdz_start).count();
        auto MP_schedule_mm_sdz_time_msec = chrono::duration_cast<chrono::milliseconds>(MP_schedule_mm_sdz_end - MP_schedule_mm_sdz_start).count();

        /*計算結果を各種配列に保存する*/

        /*時間測定結果を各種配列に保存する*/
        time_mm_sdz.push_back(mm_sdz_time_msec);
        time_MP_sdz.push_back(MP_mm_sdz_time_msec);
        time_MP_schedule_sdz.push_back(MP_schedule_mm_sdz_time_msec);
    }

    // 計算結果をファイル出力する
    std::ofstream ofs_calc(PATH_result_of_calc);
    // fileをopenできるか確認する
    if (!ofs_calc)
    {
        printf("@test1_mm_dnrm2 error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            ofs_calc << "------------------------------------------------------------------------------\n";
            ofs_calc << "@M[1][3]\n";
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << M2_13[i] << "," << M3_13[i] << "," << M4_13[i] << std::endl;
            }
            ofs_calc << "------------------------------------------------------------------------------\n";

            ofs_calc << "------------------------------------------------------------------------------\n";
            ofs_calc << "@M[10][10]\n";
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << M2_tt[i] << "," << M3_tt[i] << "," << M4_tt[i] << std::endl;
            }
            ofs_calc << "------------------------------------------------------------------------------\n";
        }
        else // txt出力
        {
            ofs_calc << "------------------------------------------------------------------------------\n";
            ofs_calc << "@M[1][3]\n";
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << M2_13[i] << " " << M3_13[i] << " " << M4_13[i] << std::endl;
            }
            ofs_calc << "------------------------------------------------------------------------------\n";

            ofs_calc << "------------------------------------------------------------------------------\n";
            ofs_calc << "@M[10][10]\n";
            for (int i = 0; i < trial_num; i++)
            {
                ofs_calc << setprecision(16) << M2_tt[i] << " " << M3_tt[i] << " " << M4_tt[i] << std::endl;
            }
            ofs_calc << "------------------------------------------------------------------------------\n";
        }
    }
    ofs_calc.close();

    // 時間測定結果をファイル出力する
    std::ofstream ofs_time(PATH_result_of_time);
    // fileをopenできるか確認する
    if (!ofs_time)
    {
        printf("@test1_mm_dnrm2 error:: \"%s\" could not open.", PATH_result_of_time.c_str());
    }
    else
    {
        if (mode == 'c') // csv出力
        {
            ofs_time << "Normal"
                     << ","
                     << "OpenMP"
                     << "OpenMP with schedule" << endl;
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_sdz[i] << "," << time_MP_sdz[i] << "," << time_MP_schedule_sdz[i] << std::endl;
            }
        }
        else // txt出力
        {
            for (int i = 0; i < trial_num; i++)
            {
                ofs_time << time_mm_sdz[i] << " " << time_MP_sdz[i] << " " << time_MP_schedule_sdz[i] << std::endl;
            }
        }
    }
    ofs_time.close();

    /*---------------------- -メモリの開放-----------------------*/
    for (int i = 0; i < mat_dim; i++)
        delete[] M1[i];
    delete[] M1;

    for (int i = 0; i < mat_dim; i++)
        delete[] M2[i];
    delete[] M2;

    for (int i = 0; i < mat_dim; i++)
        delete[] M3[i];
    delete[] M3;
};

void test_isoA_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 初期状態行列の用意
    // dim_A, dim_Bは4とする
    dim_A = 4;
    dim_B = 4;

    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.1);
    uniform_real_distribution<> rand2(0, 0.1);

    double **M1 = new double *[dim_A];
    for (int i = 0; i < dim_A; i++)
        M1[i] = new double[dim_A];

    double **M2 = new double *[dim_B];
    for (int i = 0; i < dim_B; i++)
        M2[i] = new double[dim_B];

    double **M3 = new double *[dim_B];
    for (int i = 0; i < dim_B; i++)
        M3[i] = new double[dim_B];

    // M2とM3を初期化するための行列
    double **M4 = new double *[dim_B];
    for (int i = 0; i < dim_B; i++)
        M4[i] = new double[dim_B];

    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            M1[i][j] = rand1(mt);
            M2[i][j] = rand2(mt);
            M3[i][j] = M2[i][j]; // testのためにM3とM2は等価な行列とする
            M4[i][j] = M2[i][j];
        }
    }

    // 疎行列ｰ密行列積の実施と計算時間の測定を行う
    chrono::system_clock::time_point mmprod_start, MP_mmprod_start, mmprod_end, MP_mmprod_end;
    cout << "\n\n Test : isoA_mmprod\n";
    cout << "====================================================\n";
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mmprod;
    vector<double> time_MP_mmprod;

    // trial_num回だけテストを実施する
    for (int i = 0; i < trial_num; i++)
    {
        // About isoA_mmprod
        mmprod_start = chrono::system_clock::now();
        isoA_mmprod(dim_A, dim_B, M1, M2);
        mmprod_end = chrono::system_clock::now();

        // About MP_isoA_mmprod
        MP_mmprod_start = chrono::system_clock::now();
        MP_isoA_mmprod(dim_A, dim_B, M1, M3);
        MP_mmprod_end = chrono::system_clock::now();

        // 時間計測結果をミリ秒に変換
        auto mmprod_time_msec = chrono::duration_cast<chrono::milliseconds>(mmprod_end - mmprod_start).count();
        auto MP_mmprod_msec = chrono::duration_cast<chrono::milliseconds>(MP_mmprod_end - MP_mmprod_start).count();

        // 時間計測結果を各種配列に保存する
        time_mmprod.push_back(mmprod_time_msec);
        time_MP_mmprod.push_back(MP_mmprod_msec);

        // 計算結果をファイル出力する
        std::ofstream ofs_calc(PATH_result_of_calc, std::ios::app);
        // fileをopenできる確認する
        if (!ofs_calc)
        {
            printf("@test_isoA_mmprod error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
        }
        else
        {
            if (mode == 'c') // csv出力
            {
                ofs_calc << "@isoA_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < dim_A; n++)
                {
                    for (int m = 0; m < dim_B; m++)
                    {
                        ofs_calc << setprecision(16) << M2[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";

                ofs_calc << "@MP_isoA_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < dim_A; n++)
                {
                    for (int m = 0; m < dim_B; m++)
                    {
                        ofs_calc << setprecision(16) << M3[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";
            }
        }
        ofs_calc.close();

        // 行列M2,M3の初期化
        for (int n = 0; n < dim_A; n++)
        {
            for (int m = 0; m < dim_B; m++)
            {
                M2[n][m] = M4[n][m];
                M3[n][m] = M4[n][m];
            }
        }
    }
    // relase memory
    for (int i = 0; i < dim_A; i++)
        delete[] M1[i];
    delete[] M1;

    for (int i = 0; i < dim_B; i++)
        delete[] M2[i];
    delete[] M2;

    for (int i = 0; i < dim_B; i++)
        delete[] M3[i];
    delete[] M3;

    for (int i = 0; i < dim_B; i++)
        delete[] M4[i];
    delete[] M4;
}

void test_isoB_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 初期状態行列の用意
    // dim_A, dim_Bは4とする
    dim_A = 4;
    dim_B = 4;

    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.1);
    uniform_real_distribution<> rand2(0, 0.1);

    double **M1 = new double *[dim_A];
    for (int i = 0; i < dim_A; i++)
        M1[i] = new double[dim_A];

    double **M2 = new double *[dim_B];
    for (int i = 0; i < dim_B; i++)
        M2[i] = new double[dim_B];

    double **M3 = new double *[dim_B];
    for (int i = 0; i < dim_B; i++)
        M3[i] = new double[dim_B];

    // M2とM3を初期化するための行列
    double **M4 = new double *[dim_B];
    for (int i = 0; i < dim_B; i++)
        M4[i] = new double[dim_B];

    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            M1[i][j] = rand1(mt);
            M2[i][j] = rand2(mt);
            M3[i][j] = M2[i][j]; // testのためにM3とM2は等価な行列とする
            M4[i][j] = M2[i][j];
        }
    }

    // 疎行列ｰ密行列積の実施と計算時間の測定を行う
    chrono::system_clock::time_point mmprod_start, MP_mmprod_start, mmprod_end, MP_mmprod_end;
    cout << "\n\n Test : isoB_mmprod\n";
    cout << "====================================================\n";
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mmprod;
    vector<double> time_MP_mmprod;

    // trial_num回だけテストを実施する
    for (int i = 0; i < trial_num; i++)
    {
        // About isoA_mmprod
        mmprod_start = chrono::system_clock::now();
        isoB_mmprod(dim_A, dim_B, M1, M2);
        mmprod_end = chrono::system_clock::now();

        // About MP_isoA_mmprod
        MP_mmprod_start = chrono::system_clock::now();
        MP_isoB_mmprod(dim_A, dim_B, M1, M3);
        MP_mmprod_end = chrono::system_clock::now();

        // 時間計測結果をミリ秒に変換
        auto mmprod_time_msec = chrono::duration_cast<chrono::milliseconds>(mmprod_end - mmprod_start).count();
        auto MP_mmprod_msec = chrono::duration_cast<chrono::milliseconds>(MP_mmprod_end - MP_mmprod_start).count();

        // 時間計測結果を各種配列に保存する
        time_mmprod.push_back(mmprod_time_msec);
        time_MP_mmprod.push_back(MP_mmprod_msec);

        // 計算結果をファイル出力する
        std::ofstream ofs_calc(PATH_result_of_calc, std::ios::app);
        // fileをopenできる確認する
        if (!ofs_calc)
        {
            printf("@test_isoB_mmprod error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
        }
        else
        {
            if (mode == 'c') // csv出力
            {
                ofs_calc << "@isoB_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < dim_A; n++)
                {
                    for (int m = 0; m < dim_B; m++)
                    {
                        ofs_calc << setprecision(16) << M2[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";

                ofs_calc << "@MP_isoB_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < dim_A; n++)
                {
                    for (int m = 0; m < dim_B; m++)
                    {
                        ofs_calc << setprecision(16) << M3[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";
            }
        }
        ofs_calc.close();

        // 行列M2,M3の初期化
        for (int n = 0; n < dim_A; n++)
        {
            for (int m = 0; m < dim_B; m++)
            {
                M2[n][m] = M4[n][m];
                M3[n][m] = M4[n][m];
            }
        }
    }
    // relase memory
    for (int i = 0; i < dim_A; i++)
        delete[] M1[i];
    delete[] M1;

    for (int i = 0; i < dim_B; i++)
        delete[] M2[i];
    delete[] M2;

    for (int i = 0; i < dim_B; i++)
        delete[] M3[i];
    delete[] M3;

    for (int i = 0; i < dim_B; i++)
        delete[] M4[i];
    delete[] M4;
}

void test_int_rise_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 行列計算の例として8site totalS^z = 1での(S_A^z , S_B^z) = (0,1)での計算を考える
    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.1);

    // V1は6行4列の状態行列
    double **V1 = new double *[6];
    for (int i = 0; i < 6; i++)
        V1[i] = new double[4];

    // V2は4行6列の状態行列
    double **V2 = new double *[4];
    for (int i = 0; i < 4; i++)
        V2[i] = new double[6];

    // V3は4行6列の状態行列
    double **V3 = new double *[4];
    for (int i = 0; i < 4; i++)
        V3[i] = new double[6];

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            V1[i][j] = rand1(mt);
        }
    }

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            V2[i][j] = 0.0;
            V3[i][j] = 0.0;
        }
    }

    // 昇降演算子の情報を記録するための配列のメモリ確保&要素代入
    int *prow_ind_A = new int[3];
    int *pcol_ind_A = new int[3];
    int *prow_ind_B = new int[3];
    int *pcol_ind_B = new int[3];

    prow_ind_A[0] = 0;
    prow_ind_A[1] = 1;
    prow_ind_A[2] = 2;
    pcol_ind_A[0] = 2;
    pcol_ind_A[1] = 4;
    pcol_ind_A[2] = 5;

    prow_ind_B[0] = 1;
    prow_ind_B[1] = 2;
    prow_ind_B[2] = 3;
    pcol_ind_B[0] = 0;
    pcol_ind_B[1] = 1;
    pcol_ind_B[2] = 2;

    // 疎行列ｰ密行列積の実施と計算時間の測定を行う
    chrono::system_clock::time_point mmprod_start, MP_mmprod_start, mmprod_end, MP_mmprod_end;
    cout << "\n\n Test : int_rise_mmprod\n";
    cout << "====================================================\n";
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mmprod;
    vector<double> time_MP_mmprod;

    // trial_num回だけテストを実施する
    std::ofstream ofs_calc(PATH_result_of_calc);
    for (int i = 0; i < trial_num; i++)
    {
        // About isoA_mmprod
        mmprod_start = chrono::system_clock::now();
        int_rise_mmprod(3, 3, prow_ind_A, pcol_ind_A, prow_ind_B, pcol_ind_B, V1, V2);
        mmprod_end = chrono::system_clock::now();

        // About MP_isoA_mmprod
        MP_mmprod_start = chrono::system_clock::now();
        MP_int_rise_mmprod(3, 3, prow_ind_A, pcol_ind_A, prow_ind_B, pcol_ind_B, V1, V3);
        MP_mmprod_end = chrono::system_clock::now();

        // 時間計測結果をミリ秒に変換
        auto mmprod_time_msec = chrono::duration_cast<chrono::milliseconds>(mmprod_end - mmprod_start).count();
        auto MP_mmprod_msec = chrono::duration_cast<chrono::milliseconds>(MP_mmprod_end - MP_mmprod_start).count();

        // 時間計測結果を各種配列に保存する
        time_mmprod.push_back(mmprod_time_msec);
        time_MP_mmprod.push_back(MP_mmprod_msec);

        // 計算結果をファイル出力する
        // fileをopenできる確認する
        if (!ofs_calc)
        {
            printf("@test_int_rise_mmprod error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
        }
        else
        {
            if (mode == 'c') // csv出力
            {
                ofs_calc << "@int_rise_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < 4; n++)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        ofs_calc << setprecision(16) << V2[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";

                ofs_calc << "@MP_int_rise_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < 4; n++)
                {
                    for (int m = 0; m < 6; m++)
                    {
                        ofs_calc << setprecision(16) << V3[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";
            }
        }

        // 行列M2,M3の初期化
        for (int n = 0; n < 4; n++)
        {
            for (int m = 0; m < 6; m++)
            {
                V2[n][m] = 0.0;
                V3[n][m] = 0.0;
            }
        }
    }
    ofs_calc.close();
    // relase memory
    for (int i = 0; i < 6; i++)
        delete[] V1[i];
    delete[] V1;

    for (int i = 0; i < 4; i++)
        delete[] V2[i];
    delete[] V2;

    for (int i = 0; i < 4; i++)
        delete[] V3[i];
    delete[] V3;

    delete[] prow_ind_A;
    delete[] pcol_ind_A;
    delete[] prow_ind_B;
    delete[] pcol_ind_B;
};

void test_int_dsmn_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 行列計算の例として8site totalS^z = 1での(S_A^z , S_B^z) = (0,1)での計算を考える
    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.1);

    // V1は4行6列の状態行列
    double **V1 = new double *[4];
    for (int i = 0; i < 4; i++)
        V1[i] = new double[6];

    // V2は6行4列の状態行列
    double **V2 = new double *[6];
    for (int i = 0; i < 6; i++)
        V2[i] = new double[4];

    // V3は6行4列の状態行列
    double **V3 = new double *[6];
    for (int i = 0; i < 6; i++)
        V3[i] = new double[4];

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            V1[i][j] = rand1(mt);
        }
    }

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            V2[i][j] = 0.0;
            V3[i][j] = 0.0;
        }
    }

    // 昇降演算子の情報を記録するための配列のメモリ確保&要素代入
    int *mrow_ind_A = new int[3];
    int *mcol_ind_A = new int[3];
    int *mrow_ind_B = new int[3];
    int *mcol_ind_B = new int[3];

    mrow_ind_A[0] = 2;
    mrow_ind_A[1] = 4;
    mrow_ind_A[2] = 5;
    mcol_ind_A[0] = 0;
    mcol_ind_A[1] = 1;
    mcol_ind_A[2] = 2;

    mrow_ind_B[0] = 0;
    mrow_ind_B[1] = 1;
    mrow_ind_B[2] = 2;
    mcol_ind_B[0] = 1;
    mcol_ind_B[1] = 2;
    mcol_ind_B[2] = 3;

    // 疎行列ｰ密行列積の実施と計算時間の測定を行う
    chrono::system_clock::time_point mmprod_start, MP_mmprod_start, mmprod_end, MP_mmprod_end;
    cout << "\n\n Test : int_dsmn_mmprod\n";
    cout << "====================================================\n";
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mmprod;
    vector<double> time_MP_mmprod;

    // trial_num回だけテストを実施する
    std::ofstream ofs_calc(PATH_result_of_calc);
    for (int i = 0; i < trial_num; i++)
    {
        // About isoA_mmprod
        mmprod_start = chrono::system_clock::now();
        int_dsmn_mmprod(3, 3, mrow_ind_A, mcol_ind_A, mrow_ind_B, mcol_ind_B, V1, V2);
        mmprod_end = chrono::system_clock::now();

        // About MP_isoA_mmprod
        MP_mmprod_start = chrono::system_clock::now();
        MP_int_dsmn_mmprod(3, 3, mrow_ind_A, mcol_ind_A, mrow_ind_B, mcol_ind_B, V1, V3);
        MP_mmprod_end = chrono::system_clock::now();

        // 時間計測結果をミリ秒に変換
        auto mmprod_time_msec = chrono::duration_cast<chrono::milliseconds>(mmprod_end - mmprod_start).count();
        auto MP_mmprod_msec = chrono::duration_cast<chrono::milliseconds>(MP_mmprod_end - MP_mmprod_start).count();

        // 時間計測結果を各種配列に保存する
        time_mmprod.push_back(mmprod_time_msec);
        time_MP_mmprod.push_back(MP_mmprod_msec);

        // 計算結果をファイル出力する
        // fileをopenできる確認する
        if (!ofs_calc)
        {
            printf("@test_int_dsmn_mmprod error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
        }
        else
        {
            if (mode == 'c') // csv出力
            {
                ofs_calc << "@int_dsmn_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < 6; n++)
                {
                    for (int m = 0; m < 4; m++)
                    {
                        ofs_calc << setprecision(16) << V2[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";

                ofs_calc << "@MP_int_dsmn_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < 6; n++)
                {
                    for (int m = 0; m < 4; m++)
                    {
                        ofs_calc << setprecision(16) << V3[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";
            }
        }

        // 行列M2,M3の初期化
        for (int n = 0; n < 6; n++)
        {
            for (int m = 0; m < 4; m++)
            {
                V2[n][m] = 0.0;
                V3[n][m] = 0.0;
            }
        }
    }
    ofs_calc.close();
    // relase memory
    for (int i = 0; i < 4; i++)
        delete[] V1[i];
    delete[] V1;

    for (int i = 0; i < 6; i++)
        delete[] V2[i];
    delete[] V2;

    for (int i = 0; i < 6; i++)
        delete[] V3[i];
    delete[] V3;

    delete[] mrow_ind_A;
    delete[] mcol_ind_A;
    delete[] mrow_ind_B;
    delete[] mcol_ind_B;
};

void test_int_zz_mmprod(int trial_num, int dim_A, int dim_B, std::string PATH_result_of_calc, std::string PATH_result_of_time, char mode)
{
    // 行列計算の例として8site totalS^z = 1での(S_A^z , S_B^z) = (0,1)での計算を考える
    random_device rand;
    mt19937 mt(rand());
    uniform_real_distribution<> rand1(0, 0.1);

    // V1は6行4列の状態行列
    double **V1 = new double *[6];
    for (int i = 0; i < 6; i++)
        V1[i] = new double[4];

    // V2は6行4列の状態行列
    double **V2 = new double *[6];
    for (int i = 0; i < 6; i++)
        V2[i] = new double[4];

    // V3は6行4列の状態行列
    double **V3 = new double *[6];
    for (int i = 0; i < 6; i++)
        V3[i] = new double[4];

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            V1[i][j] = rand1(mt);
            V2[i][j] = 0.0;
            V3[i][j] = 0.0;
        }
    }

    // 昇降演算子の情報を記録するための配列のメモリ確保&要素代入
    int *sz_A = new int[6];
    int *sz_B = new int[4];

    for (int i = 0; i < 6; i++)
        sz_A[i] = 1;

    for (int i = 0; i < 4; i++)
        sz_B[i] = 1;

    // 疎行列ｰ密行列積の実施と計算時間の測定を行う
    chrono::system_clock::time_point mmprod_start, MP_mmprod_start, mmprod_end, MP_mmprod_end;
    cout << "\n\n Test : int_zz_mmprod\n";
    cout << "====================================================\n";
    // 各種内積計算の時間測定結果を格納する
    vector<double> time_mmprod;
    vector<double> time_MP_mmprod;

    // trial_num回だけテストを実施する
    std::ofstream ofs_calc(PATH_result_of_calc);
    for (int i = 0; i < trial_num; i++)
    {
        // About isoA_mmprod
        mmprod_start = chrono::system_clock::now();
        int_zz_mmprod(6, 4, sz_A, sz_B, V1, V2);
        mmprod_end = chrono::system_clock::now();

        // About MP_isoA_mmprod
        MP_mmprod_start = chrono::system_clock::now();
        MP_int_zz_mmprod(6, 4, sz_A, sz_B, V1, V3);
        MP_mmprod_end = chrono::system_clock::now();

        // 時間計測結果をミリ秒に変換
        auto mmprod_time_msec = chrono::duration_cast<chrono::milliseconds>(mmprod_end - mmprod_start).count();
        auto MP_mmprod_msec = chrono::duration_cast<chrono::milliseconds>(MP_mmprod_end - MP_mmprod_start).count();

        // 時間計測結果を各種配列に保存する
        time_mmprod.push_back(mmprod_time_msec);
        time_MP_mmprod.push_back(MP_mmprod_msec);

        // 計算結果をファイル出力する
        // fileをopenできる確認する
        if (!ofs_calc)
        {
            printf("@test_int_zz_mmprod error:: \"%s\" could not open.", PATH_result_of_calc.c_str());
        }
        else
        {
            if (mode == 'c') // csv出力
            {
                ofs_calc << "@int_zz_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < 6; n++)
                {
                    for (int m = 0; m < 4; m++)
                    {
                        ofs_calc << setprecision(16) << V2[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";

                ofs_calc << "@MP_int_zz_mmprod's result @trial_num = " << i << endl;
                ofs_calc << "--------------------------------------\n";
                for (int n = 0; n < 6; n++)
                {
                    for (int m = 0; m < 4; m++)
                    {
                        ofs_calc << setprecision(16) << V3[n][m] << " , ";
                    }
                    ofs_calc << "\n";
                }
                ofs_calc << "--------------------------------------\n";
            }
        }

        // 行列M2,M3の初期化
        for (int n = 0; n < 6; n++)
        {
            for (int m = 0; m < 4; m++)
            {
                V2[n][m] = 0.0;
                V3[n][m] = 0.0;
            }
        }
    }
    ofs_calc.close();
    // relase memory
    for (int i = 0; i < 6; i++)
        delete[] V1[i];
    delete[] V1;

    for (int i = 0; i < 6; i++)
        delete[] V2[i];
    delete[] V2;

    for (int i = 0; i < 6; i++)
        delete[] V3[i];
    delete[] V3;

    delete[] sz_A;
    delete[] sz_B;
};
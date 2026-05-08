#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>


constexpr double pi = 3.14159265358979323846;


struct wir {
    double x = 0.0;
    double y = 0.0;
    double epsilon = 0.0;
    double gamma = 0.0;

    bool operator!=(const wir& other) const {
        return x != other.x
            || y != other.y
            || epsilon != other.epsilon
            || gamma != other.gamma;
    }
};


struct ZakresProcesu {
    int start = 0;
    int end = 0;
};


double cyrkulacja_cala(double u1, double u2, double L) {
    return u1 * L - u2 * L;
}


int len(const std::vector<wir>& arr) {
    return static_cast<int>(arr.size());
}


ZakresProcesu zakres_dla_procesu(int N, int rank, int size) {
    int base = N / size;
    int rest = N % size;

    int start = rank * base + std::min(rank, rest);
    int local_count = base + (rank < rest ? 1 : 0);

    ZakresProcesu zakres;
    zakres.start = start;
    zakres.end = start + local_count;

    return zakres;
}


void przygotuj_mpi_layout(int N, int size, std::vector<int>& counts, std::vector<int>& displs) {
    counts.resize(size);
    displs.resize(size);

    for (int r = 0; r < size; ++r) {
        ZakresProcesu zakres = zakres_dla_procesu(N, r, size);

        counts[r] = static_cast<int>((zakres.end - zakres.start) * sizeof(wir));
        displs[r] = static_cast<int>(zakres.start * sizeof(wir));
    }
}


void synchronizuj_wiry(std::vector<wir>& wiry) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = len(wiry);

    std::vector<int> counts;
    std::vector<int> displs;

    przygotuj_mpi_layout(N, size, counts, displs);

    MPI_Allgatherv(
        MPI_IN_PLACE,
        0,
        MPI_DATATYPE_NULL,
        wiry.data(),
        counts.data(),
        displs.data(),
        MPI_BYTE,
        MPI_COMM_WORLD
    );
}


// redyskretyzacja
void modifyTablicaWirow(std::vector<wir>& wiry, double l_krytyczne_kwadrat) {
    int i = 1;

    while (i < len(wiry)) {
        double dx = wiry[i].x - wiry[i - 1].x;
        double dy = wiry[i].y - wiry[i - 1].y;

        double distance_squared = dx * dx + dy * dy;

        if (distance_squared > l_krytyczne_kwadrat) {
            wir nowy_element;

            nowy_element.x = 0.5 * (wiry[i].x + wiry[i - 1].x);
            nowy_element.y = 0.5 * (wiry[i].y + wiry[i - 1].y);
            nowy_element.epsilon = 0.5 * (wiry[i].epsilon + wiry[i - 1].epsilon);
            nowy_element.gamma = 0.5 * (wiry[i].gamma + wiry[i - 1].gamma);

            wiry.insert(wiry.begin() + i, nowy_element);
        } else {
            ++i;
        }
    }
}


// predkosc x
double vx(const wir& wir_1, const wir& wir_2, double delta) {
    double dy = wir_1.y - wir_2.y;
    double dx = wir_1.x - wir_2.x;

    double licznik = -0.5 * wir_2.gamma * std::sinh(2.0 * pi * dy);

    double mianownik =
        std::cosh(2.0 * pi * dy)
        - std::cos(2.0 * pi * dx)
        + delta * delta;

    return licznik / mianownik;
}


// predkosc y
double vy(const wir& wir_1, const wir& wir_2, double delta) {
    double dy = wir_1.y - wir_2.y;
    double dx = wir_1.x - wir_2.x;

    double licznik = 0.5 * wir_2.gamma * std::sin(2.0 * pi * dx);

    double mianownik =
        std::cosh(2.0 * pi * dy)
        - std::cos(2.0 * pi * dx)
        + delta * delta;

    return licznik / mianownik;
}


double vx_wypadkowa(const wir& wir_0, const std::vector<wir>& wiry, double delta) {
    double result = 0.0;
    int N = len(wiry);

    for (int i = 1; i < N; ++i) {
        double d_epsilon = wiry[i].epsilon - wiry[i - 1].epsilon;
        result += vx(wir_0, wiry[i], delta) * d_epsilon;
    }

    return result;
}


double vy_wypadkowa(const wir& wir_0, const std::vector<wir>& wiry, double delta) {
    double result = 0.0;
    int N = len(wiry);

    for (int i = 1; i < N; ++i) {
        double d_epsilon = wiry[i].epsilon - wiry[i - 1].epsilon;
        result += 0.5 * d_epsilon * vy(wir_0, wiry[i], delta);
    }

    return result;
}


void euler(std::vector<wir>& tablica_wirow, double delta, double delta_t) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = len(tablica_wirow);
    ZakresProcesu zakres = zakres_dla_procesu(N, rank, size);

    std::vector<wir> predkosci_wirow(N);

    for (int i = zakres.start; i < zakres.end; ++i) {
        predkosci_wirow[i].x = vx_wypadkowa(tablica_wirow[i], tablica_wirow, delta);
        predkosci_wirow[i].y = vy_wypadkowa(tablica_wirow[i], tablica_wirow, delta);
    }

    for (int i = zakres.start; i < zakres.end; ++i) {
        tablica_wirow[i].x += predkosci_wirow[i].x * delta_t;
        tablica_wirow[i].y += predkosci_wirow[i].y * delta_t;
    }

    synchronizuj_wiry(tablica_wirow);
}


std::vector<wir> runge_kutta(
    const std::vector<wir>& tablica_wirow,
    double delta,
    double delta_t,
    double& commTime,
    double& compTime
) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = len(tablica_wirow);
    ZakresProcesu zakres = zakres_dla_procesu(N, rank, size);

    std::vector<double> k1_x(N, 0.0);
    std::vector<double> k1_y(N, 0.0);
    std::vector<double> k2_x(N, 0.0);
    std::vector<double> k2_y(N, 0.0);
    std::vector<double> k3_x(N, 0.0);
    std::vector<double> k3_y(N, 0.0);
    std::vector<double> k4_x(N, 0.0);
    std::vector<double> k4_y(N, 0.0);

    std::vector<wir> stan_1 = tablica_wirow;
    std::vector<wir> stan_2 = tablica_wirow;
    std::vector<wir> stan_3 = tablica_wirow;
    std::vector<wir> wynik = tablica_wirow;

    double start_time = 0.0;
    double end_time = 0.0;
    double comm_time = 0.0;
    double comp_time = 0.0;

    // k1
    start_time = MPI_Wtime();

    for (int i = zakres.start; i < zakres.end; ++i) {
        k1_x[i] = vx_wypadkowa(tablica_wirow[i], tablica_wirow, delta);
        k1_y[i] = vy_wypadkowa(tablica_wirow[i], tablica_wirow, delta);

        stan_1[i].x = tablica_wirow[i].x + 0.5 * delta_t * k1_x[i];
        stan_1[i].y = tablica_wirow[i].y + 0.5 * delta_t * k1_y[i];
    }

    end_time = MPI_Wtime();
    comp_time += end_time - start_time;

    start_time = MPI_Wtime();
    synchronizuj_wiry(stan_1);
    end_time = MPI_Wtime();
    comm_time += end_time - start_time;

    // k2
    start_time = MPI_Wtime();

    for (int i = zakres.start; i < zakres.end; ++i) {
        k2_x[i] = vx_wypadkowa(stan_1[i], stan_1, delta);
        k2_y[i] = vy_wypadkowa(stan_1[i], stan_1, delta);

        stan_2[i].x = tablica_wirow[i].x + 0.5 * delta_t * k2_x[i];
        stan_2[i].y = tablica_wirow[i].y + 0.5 * delta_t * k2_y[i];
    }

    end_time = MPI_Wtime();
    comp_time += end_time - start_time;

    start_time = MPI_Wtime();
    synchronizuj_wiry(stan_2);
    end_time = MPI_Wtime();
    comm_time += end_time - start_time;

    // k3
    start_time = MPI_Wtime();

    for (int i = zakres.start; i < zakres.end; ++i) {
        k3_x[i] = vx_wypadkowa(stan_2[i], stan_2, delta);
        k3_y[i] = vy_wypadkowa(stan_2[i], stan_2, delta);

        stan_3[i].x = tablica_wirow[i].x + delta_t * k3_x[i];
        stan_3[i].y = tablica_wirow[i].y + delta_t * k3_y[i];
    }

    end_time = MPI_Wtime();
    comp_time += end_time - start_time;

    start_time = MPI_Wtime();
    synchronizuj_wiry(stan_3);
    end_time = MPI_Wtime();
    comm_time += end_time - start_time;

    // k4
    start_time = MPI_Wtime();

    for (int i = zakres.start; i < zakres.end; ++i) {
        k4_x[i] = vx_wypadkowa(stan_3[i], stan_3, delta);
        k4_y[i] = vy_wypadkowa(stan_3[i], stan_3, delta);

        wynik[i].x = tablica_wirow[i].x
            + delta_t * (k1_x[i] + 2.0 * k2_x[i] + 2.0 * k3_x[i] + k4_x[i]) / 6.0;

        wynik[i].y = tablica_wirow[i].y
            + delta_t * (k1_y[i] + 2.0 * k2_y[i] + 2.0 * k3_y[i] + k4_y[i]) / 6.0;
    }

    end_time = MPI_Wtime();
    comp_time += end_time - start_time;

    start_time = MPI_Wtime();
    synchronizuj_wiry(wynik);
    end_time = MPI_Wtime();
    comm_time += end_time - start_time;

    if (rank == 0) {
        commTime += comm_time;
        compTime += comp_time;
    }

    return wynik;
}


double moment_x(const std::vector<wir>& tablica_wirow) {
    double result = 0.0;

    for (int i = 1; i < len(tablica_wirow); ++i) {
        double d_epsilon = tablica_wirow[i].epsilon - tablica_wirow[i - 1].epsilon;
        result += tablica_wirow[i].gamma * tablica_wirow[i].x * d_epsilon;
    }

    return result;
}


double moment_y(const std::vector<wir>& tablica_wirow) {
    double result = 0.0;

    for (int i = 1; i < len(tablica_wirow); ++i) {
        double d_epsilon = tablica_wirow[i].epsilon - tablica_wirow[i - 1].epsilon;
        result += tablica_wirow[i].gamma * tablica_wirow[i].y * d_epsilon;
    }

    return result;
}


double bezwladnosc(const std::vector<wir>& tablica_wirow) {
    double result = 0.0;

    for (int i = 1; i < len(tablica_wirow); ++i) {
        double d_epsilon = tablica_wirow[i].epsilon - tablica_wirow[i - 1].epsilon;

        double x = tablica_wirow[i].x;
        double y = tablica_wirow[i].y;

        result += tablica_wirow[i].gamma * (x * x + y * y) * d_epsilon;
    }

    return result;
}


std::string delta_do_nazwy(double delta) {
    std::ostringstream ss;
    ss << delta;
    return ss.str();
}


void zapisz_wiry(
    const std::vector<wir>& tablica_wirow,
    double delta,
    int krok
) {
    std::ostringstream filename_stream;
    filename_stream << "results/wyniki_" << delta_do_nazwy(delta) << "_" << krok << ".txt";

    std::ofstream plik(filename_stream.str());

    plik << std::setprecision(12);

    int N = len(tablica_wirow);

    for (int i = 0; i < N; ++i) {
        plik << tablica_wirow[i].x << " ";
    }

    for (int i = 0; i < N; ++i) {
        plik << tablica_wirow[i].y << " ";
    }

    plik << "\n";
}


void zapisz_predkosci(
    const std::vector<wir>& tablica_wirow,
    double L,
    double delta,
    int krok,
    int N0
) {
    std::ostringstream filename_stream;
    filename_stream << "results/predkosc_" << delta_do_nazwy(delta) << "_" << krok << ".txt";

    std::ofstream plik(filename_stream.str());

    plik << std::setprecision(12);

    for (int i = 0; i < N0; ++i) {
        wir pomoc;
        pomoc.x = i * L / static_cast<double>(N0 - 1);
        pomoc.y = 0.0;
        pomoc.epsilon = 0.0;
        pomoc.gamma = 0.0;

        plik << vx_wypadkowa(pomoc, tablica_wirow, delta) << " ";
    }

    for (int i = 0; i < N0; ++i) {
        wir pomoc;
        pomoc.x = i * L / static_cast<double>(N0 - 1);
        pomoc.y = 0.0;
        pomoc.epsilon = 0.0;
        pomoc.gamma = 0.0;

        plik << vy_wypadkowa(pomoc, tablica_wirow, delta) << " ";
    }

    plik << "\n";
}


void zapisz_diagnostyke(
    const std::vector<wir>& tablica_wirow,
    double delta,
    int krok
) {
    std::ostringstream filename_stream;
    filename_stream << "results/diagnostyka_" << delta_do_nazwy(delta) << ".txt";

    bool dopisz_naglowek = !std::filesystem::exists(filename_stream.str());

    std::ofstream plik(filename_stream.str(), std::ios::app);

    if (dopisz_naglowek) {
        plik << "krok liczba_wirow moment_x moment_y bezwladnosc\n";
    }

    plik << krok << " "
         << len(tablica_wirow) << " "
         << moment_x(tablica_wirow) << " "
         << moment_y(tablica_wirow) << " "
         << bezwladnosc(tablica_wirow) << "\n";
}


void warunek_poczatkowy(
    std::vector<wir>& tablica_wirow,
    double L,
    double gestoscCyrkulacji,
    double l_docelowe
) {
    int N = len(tablica_wirow);

    std::mt19937 gen(12345);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    double amplituda = l_docelowe * 0.01;

    for (int i = 0; i < N; ++i) {
        tablica_wirow[i].epsilon = static_cast<double>(i) / static_cast<double>(N - 1);

        tablica_wirow[i].x = i * L / static_cast<double>(N - 1);
        tablica_wirow[i].y = dis(gen) * amplituda;

        if (i == 0 || i == N - 1) {
            tablica_wirow[i].y = 0.0;
        }

        tablica_wirow[i].gamma = gestoscCyrkulacji;
    }
}


void glowna_petla(
    double L,
    double u1,
    double u2,
    double czas_symulacji,
    double delta
) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (delta <= 0.0 || L <= 0.0 || czas_symulacji <= 0.0) {
        if (rank == 0) {
            std::cerr << "Niepoprawne parametry symulacji.\n";
        }

        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    double l_docelowe = 0.05 * delta;
    int N = static_cast<int>(L / l_docelowe);

    if (N < 2) {
        N = 2;
    }

    int N0 = N;

    double epsilonRange = 1.0;
    double cyrkulacja = cyrkulacja_cala(u1, u2, L);

    if (std::abs(cyrkulacja) < 1e-14) {
        if (rank == 0) {
            std::cerr << "Cyrkulacja jest zbyt mala.\n";
        }

        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    double gestoscCyrkulacji = cyrkulacja / epsilonRange;
    double delta_t = 0.250 * L * l_docelowe / std::abs(cyrkulacja);

    double l_docelowe_kwadrat = l_docelowe * l_docelowe;
    double l_krytyczne_kwadrat = 1.25 * l_docelowe_kwadrat;

    int liczba_krokow = static_cast<int>(czas_symulacji / delta_t);

    if (rank == 0) {
        std::filesystem::create_directories("results");

        std::cout << "liczba procesow MPI: " << size << "\n";
        std::cout << "liczba krokow: " << liczba_krokow << "\n";

        std::cout << "delta = " << delta << "\n"
                  << "N = " << N << "\n"
                  << "G = " << cyrkulacja << "\n"
                  << "g = " << gestoscCyrkulacji << "\n"
                  << "dt = " << delta_t << "\n"
                  << "dl = " << l_docelowe << "\n"
                  << "L = " << L << "\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<wir> tablica_wirow(N);

    if (rank == 0) {
        warunek_poczatkowy(tablica_wirow, L, gestoscCyrkulacji, l_docelowe);
    }

    MPI_Bcast(
        tablica_wirow.data(),
        static_cast<int>(tablica_wirow.size() * sizeof(wir)),
        MPI_BYTE,
        0,
        MPI_COMM_WORLD
    );

    double commTime = 0.0;
    double compTime = 0.0;

    for (int krok = 0; krok < liczba_krokow; ++krok) {
        modifyTablicaWirow(tablica_wirow, l_krytyczne_kwadrat);

        if (rank == 0 && krok % 10 == 0) {
            std::cout << "krok: " << krok << "\n";
            std::cout << "    n = " << len(tablica_wirow) << "\n";
            std::cout << "    communication: " << commTime / 10.0 << " s\n";
            std::cout << "    computation:   " << compTime / 10.0 << " s\n";

            compTime = 0.0;
            commTime = 0.0;
        }

        tablica_wirow = runge_kutta(tablica_wirow, delta, delta_t, commTime, compTime);

        if (rank == 0 && krok % 10 == 0) {
            zapisz_wiry(tablica_wirow, delta, krok);
            zapisz_predkosci(tablica_wirow, L, delta, krok, N0);
            zapisz_diagnostyke(tablica_wirow, delta, krok);
        }
    }
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double L = 1.0;
    double u1 = 2.0;
    double u2 = 1.0;
    double czas_symulacji = 10.0;
    double delta = 0.1;

    if (argc == 6) {
        L = std::atof(argv[1]);
        u1 = std::atof(argv[2]);
        u2 = std::atof(argv[3]);
        czas_symulacji = std::atof(argv[4]);
        delta = std::atof(argv[5]);
    } else if (argc != 1) {
        if (rank == 0) {
            std::cerr << "Uzycie:\n";
            std::cerr << "  mpirun -np 4 ./vortex_sheet\n";
            std::cerr << "  mpirun -np 4 ./vortex_sheet L u1 u2 czas delta\n";
        }

        MPI_Finalize();
        return 1;
    }

    glowna_petla(L, u1, u2, czas_symulacji, delta);

    if (rank == 0) {
        std::cout << "done!\n";
    }

    MPI_Finalize();

    return 0;
}

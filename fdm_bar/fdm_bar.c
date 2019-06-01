/*
 * Программа предназначена для решения нестационарного уравнения теплопродоности
 * для прямоугольной пластины методом конечных разностей
 * с помощью неявной разностной схемы.
 * Производится расчет узлов сетки на временных слоях и выводится график цветовой гаммы
 * распределения температуры по пластине в каждый момент времени.
 *
 * Компиляция: gcc fdm_bar.c -lgsl -o fdm_bar
 * Запуск: ./fdm_bar
 * */

#include "fdm_bar.h"
//#include "../../development-of-software-systems-BMSTU-labs-/GaussMultithread/GaussMultithread.h"

// Установка краевых условий
void initialization(gsl_vector *v) {
    for (int i = 0; i < vector_N; i++) {
        gsl_vector_set(v, i, init_condition);
    }
    boundary_condition(v);
}

/*
 * Функция производит расчет температуры методом конечных разностей
 * с использованием неявной разностной схемы
 */
void FDM(gsl_vector *b, gsl_matrix *A) {
    int s;
    gsl_vector *x = gsl_vector_alloc(vector_N);
    gsl_permutation *p = gsl_permutation_alloc(vector_N);
    get_LAE(A, b);
    print_LAE(b, A);
//    print_temperature(b);
//    gauss_solve(A->data, b->data, x->data, vector_N, vector_N);
    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_solve(A, p, b, x);
    printf("\n");
    print_temperature(x);
    create_plot(b);
    create_plot(x);
    for (int j = 0; j < model_time; j++) {
        boundary_condition(x);
//        gauss_solve(A->data, x->data, x->data, vector_N, vector_N);
//        print_LAE(x, A);
//        print_temperature(x);
//        printf("\n");
        gsl_linalg_LU_svx(A, p, x);
        print_temperature(x);
        create_plot(x);
    }
}

// Отладочная печать СЛАУ
void print_LAE(const gsl_vector *b, const gsl_matrix *A) {
    for (int k = 0, v = 0, g = 0; k < vector_N; k++) {
        for (int i = 0; i < vector_N; i++) {
            printf("%.0f ", A->data[k * vector_N + i]);
        }
        printf("\tT[%d][%d]", v, g);
        printf("\t%.2f", b->data[k]);
        putchar('\n');
        if ((k + 1) % x_nodes == 0) {
            v++;
            g = 0;
        } else {
            g++;
        }
    }
    putchar('\n');
}

// Печать таблицы распределения температуры
void print_temperature(gsl_vector *t) {
    for (int i = y_nodes - 1; i >= 0; i--) {
        for (int j = 0; j < x_nodes; j++) {
            printf("%4.0f ", t->data[i * x_nodes + j]);
        }
        putchar('\n');
    }
}


/*
 * Составление СЛАУ с уравнениями для внутренних узлов сетки, а также
 * уравненияи для гранчных узлов стеки, на которые наложены
 * граничные условия 2 рода
 */

void get_LAE(gsl_matrix *A, gsl_vector *b) {
    float dx = 8.0 / (x_nodes - 1);
    float dy = 6.0 / (y_nodes - 1);
    for (int i = 0; i < vector_N; i++) {
        float xc = 0, yc = 0;
        for (int j = 0; j < vector_N; j++) {
            if (i == j && xc + yc < 11 + dx) {
                if (gsl_vector_get(b, i) == 0 && xc > 0.1 && yc > 0.1 && xc < 7.9 &&
                    yc < 5.9)
                {
                    gsl_matrix_set(A, i, j + 1, a / (dx * dx));
                    gsl_matrix_set(A, i, j, ((-2 * a / (dx * dx)) - (2 * a / (dy * dy)) - 1.0 / dt));
                    gsl_matrix_set(A, i, j - 1, a / (dx * dx));
                    gsl_matrix_set(A, i, j + x_nodes, a / (dy * dy));
                    gsl_matrix_set(A, i, j - x_nodes, a / (dy * dy));

                } else if (yc == 0 && gsl_vector_get(b, i) == 0) {
                    gsl_matrix_set(A, i, j, a / (dy));
                    gsl_matrix_set(A, i, i + x_nodes, -a / (dy));
                }

            }
            if ((xc + dx) >= 8.001) {
                xc = 0.0;
                yc += dy;
            } else
                xc += dx;
        }
    }
}

// Установка граничных условий 1 и 2 рода
void boundary_condition(gsl_vector *v) {

    float x = 0, y = 0;
    float dx =  8.0 / (x_nodes - 1);
    float dy =  6.0 / (y_nodes - 1);

    for (int i = 0; i < (int) x_nodes * y_nodes; i++) {
        if (y == 0) {
            gsl_vector_set(v, i, 200);
        }

        if ((abs(x - 9) < dx && y <= 3) || (abs(y - 7) < dy && x <= 5)) {
            gsl_vector_set(v, i, corner);
        }

        if (x == 0) {
            gsl_vector_set(v, i, left_bound);
        }

        if (abs(x + y - 12) < dx / 2 && x > 5 && y > 3) {
            //printf("YAY\n");
            gsl_vector_set(v, i, champher);
        }
        if (x + y > 11 + dx) {
            gsl_vector_set(v, i, 0);
        }

        x += dx;
        if (x >= 8.001) {
            x = 0.0;
            y += dy;
        }
        if (y >= 6.001) {
            x = 0;
            y = 0;
        }

    }
/*
    for (int i = 0; i < N; i++) {
        //printf("[%.3f; %.3f] %g\n", x, y, gsl_vector_get(v, i));
    }*/
    putchar('\n');
}

// Создание цветовой гаммы распределения температуры в пластине
void create_plot(gsl_vector *tepmerature) {
    float dx = (float) 8 / (x_nodes - 1);
    float dy = (float) 6 / (y_nodes - 1);
    fprintf(gp, "set cbrange [0:200]\n"
                "set yrange [0:6]\n"
                "set xrange [0:8]\n"
                "set pm3d scansforward ftriangles map interpolate 10,10\n");
    fprintf(gp, "splot '-'\n");
    for (int j = 0; j < y_nodes; ++j) {
        for (int i = 0; i < x_nodes; i++) {
            fprintf(gp, "%-15f %-15f %-15g\n", i * dx, j * dy, tepmerature->data[j * x_nodes + i]);
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");
    fprintf(gp, "pause 0.5\n");
    fflush(gp);
}

int main(int argc, char *argv[]) {
    gp = popen("gnuplot -persist", "w");
    gsl_vector *init = gsl_vector_alloc(vector_N);
    initialization(init);
    gsl_matrix *m = gsl_matrix_alloc(vector_N, vector_N);
    gsl_matrix_set_identity(m);
    FDM(init, m);
    pclose(gp);
    return 0;
}

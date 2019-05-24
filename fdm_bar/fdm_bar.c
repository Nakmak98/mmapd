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
    get_LAE(A);
    print_LAE(b, A);
    print_temperature(b);
    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_solve(A, p, b, x);
    print_temperature(x);
    create_plot(b);
    create_plot(x);
    for (int j = 0; j < model_time; j++) {
        boundary_condition(x);
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
            printf("%.0f ", t->data[i * x_nodes + j]);
        }
        putchar('\n');
    }
}


/*
 * Составление СЛАУ с уравнениями для внутренних узлов сетки, а также
 * уравненияи для гранчных узлов стеки, на которые наложены
 * граничные условия 2 рода
 */
void get_LAE(gsl_matrix *A) {
    for (int k = 1; k < y_nodes - 1; k++) {
        for (int i = 1; i < x_nodes - 1; i++) {
            gsl_matrix_set(A, i + x_nodes * k, i + 1 + x_nodes * k, -a / (dy * dy));
            gsl_matrix_set(A, i + x_nodes * k, i + x_nodes * k,
                           ((2 * a / (dx * dx)) + (2 * a / (dy * dy)) + 1.0 / dt));
            gsl_matrix_set(A, i + x_nodes * k, i - 1 + x_nodes * k, -a / (dx * dx));
            gsl_matrix_set(A, i + x_nodes * k, i + x_nodes * (k + 1), -a / (dy * dy));
            gsl_matrix_set(A, i + x_nodes * k, i + x_nodes * (k - 1), -a / (dy * dy));
            if (k == y_nodes - 2) {
                gsl_matrix_set(A, (i - 1) + x_nodes * (k + 1), (i - 1) + x_nodes * (k + 1), a / dy);
                gsl_matrix_set(A, (i - 1) + x_nodes * (k + 1), (i - 1) + x_nodes * k, -a / dy);
                gsl_matrix_set(A, (i) + x_nodes * (k + 1), (i) + x_nodes * (k + 1), a / dy);
                gsl_matrix_set(A, (i) + x_nodes * (k + 1), (i) + x_nodes * k, -a / dy);
                if (i == x_nodes - 2) {
                    gsl_matrix_set(A, (i + 1) + x_nodes * (k + 1), (i + 1) + x_nodes * (k + 1), a / dy);
                    gsl_matrix_set(A, (i + 1) + x_nodes * (k + 1), (i + 1) + x_nodes * (k), -a / dy);
                }
            }

            if ((i == x_nodes - 2) && (k >= (y_nodes) / 2)) {
                gsl_matrix_set(A, (i + 1) + x_nodes * k, (i + 1) + x_nodes * k, a / (dx / dx));
                gsl_matrix_set(A, (i + 1) + x_nodes * k, i + x_nodes * k, -a / (dx / dx));
            }
        }
    }
}
// Установка граничных условий 1 и 2 рода
void boundary_condition(gsl_vector *x) {

    for (int j = 0; j < y_nodes; j++) {
        gsl_vector_set(x, j * x_nodes, left_bound);

        if (j > 0)
            gsl_vector_set(x, j * x_nodes - 1, j <= y_nodes / 2 ?
                                               bottom_right_bound : top_right_bound);
        for (int i = 0; i < x_nodes; i++) {
            if (j == 0)
                gsl_vector_set(x, i, bottom_bound);
            if (j == y_nodes - 1)
                gsl_vector_set(x, (j) * x_nodes + i, top_bound);
        }
    }
}

// Создание цветовой гаммы распределения температуры в пластине
void create_plot(gsl_vector *tepmerature) {
    float d_x = 8.0 / (x_nodes - 1.0);
    float d_y = 3.0 / (y_nodes - 1.0);
    fprintf(gp, "set cbrange [0:400]\n"
                "set yrange [0:3]\n"
                "set xrange [0:8]\n"
                "set pm3d scansforward ftriangles map interpolate 10,10\n");
    fprintf(gp, "splot '-'\n");
    for (int j = 0; j < y_nodes; j++) {
        for (int i = 0; i < x_nodes; i++) {
            fprintf(gp, "%-15g %-15g %-15g\n", i * d_x, j * d_y, (tepmerature->data[j * x_nodes + i]));
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");
    fprintf(gp, "pause 1\n");
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
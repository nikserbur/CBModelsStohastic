#include <Rcpp.h>
#include <cmath>
#include <random>
using namespace Rcpp;

// Вспомогательная функция для вычисления ранга.
// Совпадающие значения получают средний ранг, как rank() в R (ties.method = "average"),
// иначе константный признак получал бы ранги 1..n в порядке строк.
NumericVector rank_cpp(const NumericVector& x) {
  int n = x.size();
  NumericVector ranks(n);
  std::vector<std::pair<double, int>> values;

  for (int i = 0; i < n; ++i) {
    values.push_back(std::make_pair(x[i], i));
  }

  std::sort(values.begin(), values.end());

  int i = 0;
  while (i < n) {
    int j = i;
    while (j + 1 < n && values[j + 1].first == values[i].first) ++j;
    double average_rank = (i + j) / 2.0 + 1.0;
    for (int m = i; m <= j; ++m) {
      ranks[values[m].second] = average_rank;
    }
    i = j + 1;
  }

  return ranks;
}

// Случайные нелинейные проекции копула-преобразования (Lopez-Paz, Hennig, Schölkopf, 2013):
// u = rank(x) / n, phi_i(u) = sin(s / 2 * (w_i * u + b_i)), (w_i, b_i) ~ N(0, 1).
// Веса w_i, b_i выбираются один раз на признак и общие для всех наблюдений;
// s / 2 — масштаб s / ncol([u, 1]) из оригинального кода RDC (s = 1/6).
NumericMatrix rdc_random_features(const NumericVector& ranks, int k, double s,
                                  std::mt19937& gen) {
  int n = ranks.size();
  std::normal_distribution<> norm_dist(0.0, 1.0);
  NumericMatrix features(n, k);
  for (int i = 0; i < k; ++i) {
    double w = norm_dist(gen);
    double b = norm_dist(gen);
    for (int j = 0; j < n; ++j) {
      features(j, i) = std::sin(s / 2.0 * (w * ranks[j] / n + b));
    }
  }
  return features;
}

// Наибольшая каноническая корреляция (0, если cancor() завершился ошибкой)
double max_cancor(Function& cancor_r, const NumericMatrix& x, const NumericMatrix& y) {
  try {
    List cca_result = cancor_r(x, y);
    NumericVector correlations = cca_result["cor"];
    if (correlations.size() > 0) {
      return *std::max_element(correlations.begin(), correlations.end());
    }
  } catch (...) {
  }
  return 0.0;
}

// [[Rcpp::export]]
List cancor_cpp(const NumericMatrix& x, const NumericMatrix& y) {
  int n = x.nrow();
  int p = x.ncol();
  int q = y.ncol();

  if (n <= 1 || p == 0 || q == 0) {
    return List::create(
      _["cor"] = NumericVector::create(0.0),
      _["xcoef"] = NumericMatrix(p, 1),
      _["ycoef"] = NumericMatrix(q, 1)
    );
  }

  // Центрируем данные
  NumericMatrix x_centered(n, p);
  NumericMatrix y_centered(n, q);

  for (int j = 0; j < p; ++j) {
    double mean_x = 0.0;
    for (int i = 0; i < n; ++i) {
      mean_x += x(i, j);
    }
    mean_x /= n;

    for (int i = 0; i < n; ++i) {
      x_centered(i, j) = x(i, j) - mean_x;
    }
  }

  for (int j = 0; j < q; ++j) {
    double mean_y = 0.0;
    for (int i = 0; i < n; ++i) {
      mean_y += y(i, j);
    }
    mean_y /= n;

    for (int i = 0; i < n; ++i) {
      y_centered(i, j) = y(i, j) - mean_y;
    }
  }

  // Вычисляем ковариационные матрицы
  NumericMatrix Sxx(p, p);
  NumericMatrix Syy(q, q);
  NumericMatrix Sxy(p, q);

  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += x_centered(k, i) * x_centered(k, j);
      }
      Sxx(i, j) = sum / (n - 1);
    }
  }

  for (int i = 0; i < q; ++i) {
    for (int j = 0; j < q; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += y_centered(k, i) * y_centered(k, j);
      }
      Syy(i, j) = sum / (n - 1);
    }
  }

  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < q; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += x_centered(k, i) * y_centered(k, j);
      }
      Sxy(i, j) = sum / (n - 1);
    }
  }

  // Проверка на вырожденность
  bool singular = false;
  for (int i = 0; i < p; ++i) {
    if (Sxx(i, i) < 1e-12) {
      singular = true;
      break;
    }
  }
  for (int i = 0; i < q; ++i) {
    if (Syy(i, i) < 1e-12) {
      singular = true;
      break;
    }
  }

  if (singular) {
    return List::create(
      _["cor"] = NumericVector::create(0.0),
      _["xcoef"] = NumericMatrix(p, 1),
      _["ycoef"] = NumericMatrix(q, 1)
    );
  }

  // Фолбэк на R cancor для сложных вычислений
  Function cancor_r("cancor");
  List result = cancor_r(x, y);

  return result;
}

// [[Rcpp::export]]
List RDCSelection_cpp(const DataFrame& data,
                      Nullable<String> target = R_NilValue,
                      double threshold = 0.1,
                      int k_features = 20,
                      int num_permutations = 500) {

  CharacterVector feature_names = data.names();
  int n_features = feature_names.size();
  CharacterVector features;
  NumericVector target_vector;
  bool has_target = false;

  if (target.isNotNull()) {
    String target_name = as<String>(target);
    has_target = true;

    // Находим целевую переменную
    for (int i = 0; i < n_features; ++i) {
      if (feature_names[i] == target_name) {
        target_vector = data[i];
        break;
      }
    }

    // Создаем список признаков без целевой переменной
    for (int i = 0; i < n_features; ++i) {
      if (feature_names[i] != target_name) {
        features.push_back(feature_names[i]);
      }
    }
  } else {
    features = feature_names;
  }

  int n_obs = 0;
  if (features.size() > 0) {
    NumericVector first_feature = data[as<std::string>(features[0])];
    n_obs = first_feature.size();
  }

  NumericVector rdc_scores(features.size());
  NumericVector p_values(features.size());
  CharacterVector rdc_names(features.size());

  // Масштаб случайных проекций (s = 1/6 в Lopez-Paz et al., 2013)
  const double RDC_SCALE = 1.0 / 6.0;

  // Быстрые C++ генераторы с фиксированным seed
  std::mt19937 gen(123);
  std::uniform_int_distribution<> int_dist;

  // Объявляем R функции один раз для всех циклов
  Function cancor_r("cancor");

  for (int feat_idx = 0; feat_idx < features.size(); ++feat_idx) {
    rdc_names[feat_idx] = features[feat_idx];

    NumericVector x = data[as<std::string>(features[feat_idx])];
    NumericVector y;

    if (has_target) {
      y = target_vector;
    } else {
      // Выбираем случайный другой признак (как в R: sample(features[features != feat], 1))
      CharacterVector other_features;
      for (int i = 0; i < features.size(); ++i) {
        if (i != feat_idx) {
          other_features.push_back(features[i]);
        }
      }
      if (other_features.size() > 0) {
        int_dist = std::uniform_int_distribution<>(0, other_features.size() - 1);
        int random_idx = int_dist(gen);
        y = data[as<std::string>(other_features[random_idx])];
      } else {
        y = x; // Fallback если только один признак
      }
    }

    // Обработка NA/NaN
    NumericVector x_clean(n_obs);
    NumericVector y_clean(n_obs);

    for (int i = 0; i < n_obs; ++i) {
      if (std::isnan(x[i]) || std::isinf(x[i])) {
        x_clean[i] = 0.0;
      } else {
        x_clean[i] = x[i];
      }

      if (std::isnan(y[i]) || std::isinf(y[i])) {
        y_clean[i] = 0.0;
      } else {
        y_clean[i] = y[i];
      }
    }

    // Копула-преобразование: ранги
    NumericVector x_rank = rank_cpp(x_clean);
    NumericVector y_rank = rank_cpp(y_clean);

    int k = std::min(k_features, n_obs - 1);

    // Константный признак (все ранги равны) не несёт информации о зависимости
    bool is_constant =
      *std::min_element(x_rank.begin(), x_rank.end()) == *std::max_element(x_rank.begin(), x_rank.end()) ||
      *std::min_element(y_rank.begin(), y_rank.end()) == *std::max_element(y_rank.begin(), y_rank.end());
    if (is_constant || k < 1) {
      rdc_scores[feat_idx] = 0.0;
      p_values[feat_idx] = (num_permutations > 0) ? 1.0 : NA_REAL;
      continue;
    }

    // Случайные признаки RDC с общими для всех наблюдений весами
    NumericMatrix x_matrix = rdc_random_features(x_rank, k, RDC_SCALE, gen);
    NumericMatrix y_matrix = rdc_random_features(y_rank, k, RDC_SCALE, gen);

    // RDC — наибольшая каноническая корреляция между наборами признаков
    double canonical_corr = max_cancor(cancor_r, x_matrix, y_matrix);
    rdc_scores[feat_idx] = canonical_corr;

    if (num_permutations > 0) {
      // Перестановочный тест независимости: строки признаков y переставляются
      // (те же веса, что и для наблюдаемой статистики), p = (1 + #{RDC_perm >= RDC}) / (1 + B)
      std::vector<int> row_order(n_obs);
      for (int j = 0; j < n_obs; ++j) row_order[j] = j;

      int count_greater_equal = 0;
      NumericMatrix y_matrix_perm(n_obs, k);
      for (int perm = 0; perm < num_permutations; ++perm) {
        std::shuffle(row_order.begin(), row_order.end(), gen);
        for (int i = 0; i < k; ++i) {
          for (int j = 0; j < n_obs; ++j) {
            y_matrix_perm(j, i) = y_matrix(row_order[j], i);
          }
        }
        if (max_cancor(cancor_r, x_matrix, y_matrix_perm) >= canonical_corr) {
          count_greater_equal++;
        }
      }
      p_values[feat_idx] = (1.0 + count_greater_equal) / (1.0 + num_permutations);
    } else {
      p_values[feat_idx] = NA_REAL;
    }
  }

  // Присваиваем имена
  rdc_scores.names() = rdc_names;
  p_values.names() = rdc_names;

  // Отбор признаков: с перестановками — p-value <= threshold (уровень значимости),
  // без перестановок — RDC >= threshold
  CharacterVector selected_features;
  CharacterVector removed_features;

  for (int i = 0; i < features.size(); ++i) {
    bool is_selected = (num_permutations > 0) ? (p_values[i] <= threshold)
                                              : (rdc_scores[i] >= threshold);
    if (is_selected) {
      selected_features.push_back(features[i]);
    } else {
      removed_features.push_back(features[i]);
    }
  }

  // Добавляем целевую переменную если есть
  CharacterVector final_features;
  if (has_target) {
    for (int i = 0; i < selected_features.size(); ++i) {
      final_features.push_back(selected_features[i]);
    }
    final_features.push_back(as<String>(target));
  } else {
    final_features = selected_features;
  }

  // Создаем фильтрованный DataFrame
  List filtered_data_list;
  for (int i = 0; i < final_features.size(); ++i) {
    std::string feat_name = as<std::string>(final_features[i]);
    filtered_data_list.push_back(data[feat_name]);
  }
  filtered_data_list.names() = final_features;
  DataFrame filtered_data(filtered_data_list);

  return List::create(
    _["data"] = filtered_data,
    _["rdc_scores"] = rdc_scores,
    _["p_values"] = p_values,
    _["removed_features"] = removed_features,
    _["selected_features"] = final_features,
    _["threshold"] = threshold
  );
}

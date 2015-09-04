
#ifndef __KELEMEN_MLT_H
#define __KELEMEN_MLT_H

#include <vector>
#include "xor128.h"

// A Simple and Robust Mutation Strategy for the Metropolisを参照。
// Kelemen style MLT用データ構造
// Kelemen styleではパス生成に使う乱数の空間で変異させたりする。
// その一つ一つのサンプルのデータ構造。
struct PrimarySample {
  int modify_time;
  double value;
  PrimarySample( RANDOM::xor128& xor ) {
    modify_time = 0;
    value = xor.rand01();
  }
};

// Kelemen MLTにおいて、パス生成に使う各種乱数はprimary spaceからもってくる。
// PrimarySample()を通常のrand01()の代わりに使ってパス生成する。今回は普通のパストレースを使った。（双方向パストレ等も使える）
// Metropolis法なので現在の状態から次の状態へと遷移をするがこの遷移の空間が従来のワールド空間ではなく
// パスを生成するのに使われる乱数の空間になっている。
// 乱数空間における変異（Mutate()）の結果を使って再び同じようにパストレでパスを生成するとそのパスは自然に変異後のパスになっている。
struct KelemenMLT {
private:

  // LuxRenderから拝借してきた変異関数
  inline double Mutate(const double  x) {
    const double r = xor128_.rand01();
    const double s1 = 1.0 / 512.0, s2 = 1.0 / 16.0;
    const double dx = s1 / (s1 / s2 + fabsf(2.0 * r - 1.0)) - s1 / (s1 / s2 + 1.0);
    if (r < 0.5) {
      const double x1 = x + dx;
      return (x1 < 1.0) ? x1 : x1 - 1.0;
    } else {
      const double x1 = x - dx;
      return (x1 < 0.0) ? x1 + 1. : x1;
    }
  }

public:

  RANDOM::xor128 xor128_;
  
  // 論文で使われているものとほぼ同じ
  int global_time;
  int large_step;
  int large_step_time;
  int used_rand_coords;

  std::vector<PrimarySample> primary_samples;
  std::stack <PrimarySample> primary_samples_stack;

  KelemenMLT() {
    global_time = large_step = large_step_time = used_rand_coords = 0;
    //primary_samples.resize(128);
    for(int i=0;i<128;i++)
      primary_samples.push_back( PrimarySample(xor128_) );
  }
  void InitUsedRandCoords() {
    used_rand_coords = 0;
  }

  inline double rand01( void ){ return xor128_.rand01(); }
  // 通常の乱数のかわりにこれを使う
  // 論文にのっているコードとほぼ同じ
  inline double NextSample() {
    if (primary_samples.size() <= used_rand_coords) {
      for(int i=0;i<128;i++)
        primary_samples.push_back( PrimarySample(xor128_) );
      //primary_samples.resize(primary_samples.size() * 1.5); // 拡張する
    }

    if (primary_samples[used_rand_coords].modify_time < global_time) {
      if (large_step > 0) {
        primary_samples_stack.push(primary_samples[used_rand_coords]);
        primary_samples[used_rand_coords].modify_time = global_time;
        primary_samples[used_rand_coords].value = xor128_.rand01();
      } else {
        if (primary_samples[used_rand_coords].modify_time < large_step_time) {
          primary_samples[used_rand_coords].modify_time = large_step_time;
          primary_samples[used_rand_coords].value = xor128_.rand01();
        }

        while (primary_samples[used_rand_coords].modify_time < global_time - 1) {
          primary_samples[used_rand_coords].value = Mutate(primary_samples[used_rand_coords].value);
          primary_samples[used_rand_coords].modify_time ++;
        }
        primary_samples_stack.push(primary_samples[used_rand_coords]);
        primary_samples[used_rand_coords].value = Mutate(primary_samples[used_rand_coords].value);
        primary_samples[used_rand_coords].modify_time = global_time;
      }
    }

    used_rand_coords ++;
    return primary_samples[used_rand_coords - 1].value;
  }
};

#endif

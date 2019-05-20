# KalmanFilter_Example

Kalman.C

平面上单个带电粒子在垂直平面的磁场下运动

噪声来源：初始位置；$t_k$至$t_{k+1}$时的状态随机改变（摩擦，碰撞）；磁场不均匀性

能测量到每个时刻$t_k$的粒子位置信息

Kalman_1DMeasure.C

能测量到每个时刻$t_k$的粒子位置的一维信息（横坐标）





Kalman.C :X(k/k - 1)(预测)，并给粒子实际运动时的磁场加了
随机误差
Kalman_1DMeasure.C : 只测量x 坐标，进行对比
Kalman_filter.C：X(k/N)(光滑)，与X(k/k - 1)(预测) 对比
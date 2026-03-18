这是我写过的FSM代码里最简洁最工整的一个版本，实现了三维非等间距网格因式分解程函方程的高精度走时计算。

This is the most concise and well-organized version of the FSM code I've written, which achieves high-precision travel time calculations for the factorized eikonal equation on three-dimensional non-uniform grids.

你可以直接用g++编译，使用./运行代码。

You can compile it directly using g++ and run the code with ./ .

![image](https://github.com/zhangjm-geo/Fast-Sweeping-Method-FSM-for-Factorized-Eikonal-equations/blob/main/FSM_VS_FactoredFSM1.jpg)

这是相同色标下FSM 和 Factored FSM的误差对比（这么看Factored FSM貌似没有误差）。

This is a comparison of errors between FSM and Factored FSM under the same color scale (It seems that Factored FSM has no errors).

![image](https://github.com/zhangjm-geo/Fast-Sweeping-Method-FSM-for-Factorized-Eikonal-equations/blob/main/FSM_VS_FactoredFSM2.jpg)

这是FSM 和 Factored FSM各自的真实误差（Factored FSM实际还是有误差的，只是很小）。

This is the true error of FSM and Factored FSM respectively (Factored FSM actually still has errors, but they are very small).


如果遇到任何问题，欢迎通过邮件zhangjianming@ouc.edu.cn联系我。

If you encounter any issues, please feel free to contact me via email at zhangjianming@ouc.edu.cn.

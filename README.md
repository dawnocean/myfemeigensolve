# EigenSolve — 电磁波导模式求解器

基于 MFEM 的全矢量有限元波导模式求解器，支持开放型介质波导（如硅基光波导）的截面模式分析。

## 功能特性

- **全矢量公式**：采用 Nedelec 边元素（E_t）+ H1 节点元素（E_z）混合有限元空间，正确处理横向-纵向场分量耦合
- **灵活边界条件**：PEC、PMC、PML（完美匹配层）、对称/反对称边界，逐边独立配置
- **自动网格生成**：Cartesian 2D 网格，芯层边界自适应加密，支持三角形和四边形单元
- **灵活配置**：JSON 配置文件指定波长、计算域、波导几何、材料折射率、求解参数
- **场分布输出**：ParaView VTU 格式 + PNG 图片，可直接可视化模式电场分布
- **三种特征值求解器**：内置 Block Shift-and-Invert 子空间迭代，可选 SLEPc (Krylov-Schur)、ARPACK (隐式重启 Arnoldi)
- **MPI 并行**：基于 MFEM 并行框架，支持多进程运行

## 物理模型

### 波导模式特征值问题

对于 z 方向传播的波导模式 $\mathbf{E}(x,y,z) = [\mathbf{e}_t(x,y) + e_z(x,y)\hat{z}] e^{-j\beta z}$，从 Maxwell 方程出发：

$$\nabla \times \nabla \times \mathbf{E} = k_0^2 \varepsilon_r \mathbf{E}$$

在 2D 截面上离散后，消去纵向分量 $e_z$，得到关于横向场 $\mathbf{e}_t$ 的广义特征值问题：

$$\mathbf{A}\, \mathbf{e}_t = \beta^2\, \mathbf{B}_{\text{eff}}\, \mathbf{e}_t$$

其中：

| 矩阵 | 定义 | 有限元空间 |
|------|------|-----------|
| $\mathbf{K}$ | $\int (\nabla_t \times \mathbf{e}_t)(\nabla_t \times \mathbf{f}_t)\, dA$ | ND × ND (curl-curl 刚度) |
| $\mathbf{M}_\varepsilon$ | $\int \varepsilon_r\, \mathbf{e}_t \cdot \mathbf{f}_t\, dA$ | ND × ND (ε 加权质量) |
| $\mathbf{M}$ | $\int \mathbf{e}_t \cdot \mathbf{f}_t\, dA$ | ND × ND (单位质量) |
| $\mathbf{S}_{zz}$ | $\int \nabla_t e_z \cdot \nabla_t f_z\, dA$ | H1 × H1 (标量刚度) |
| $\mathbf{M}_{\varepsilon,zz}$ | $\int \varepsilon_r\, e_z f_z\, dA$ | H1 × H1 (标量质量) |
| $\mathbf{C}$ | $\int (\nabla_t \varphi) \cdot \mathbf{w}\, dA$ | ND × H1 (梯度耦合) |

复合矩阵：

$$\mathbf{A} = k_0^2 \mathbf{M}_\varepsilon - \mathbf{K}, \qquad \mathbf{D} = \mathbf{S}_{zz} - k_0^2 \mathbf{M}_{\varepsilon,zz}, \qquad \mathbf{B}_{\text{eff}} = \mathbf{M} - \mathbf{C}\mathbf{D}^{-1}\mathbf{C}^T$$

### 求解方法

支持三种特征值求解器，通过配置文件中 `solver.eigensolver` 字段选择：

#### 方法一：内置 Block Shift-and-Invert 子空间迭代（`"builtin"`，默认）

为求解靠近目标有效折射率 $n_{\text{target}}$ 的模式（对应 $\sigma = (k_0 n_{\text{target}})^2$），构建 2×2 block 移位矩阵：

$$\mathbf{S}_{\text{block}} = \begin{bmatrix} \mathbf{A} - \sigma\mathbf{M} & \mathbf{C} \\ -\sigma\mathbf{C}^T & \mathbf{D} \end{bmatrix}$$

使用 SuperLU 直接分解 $\mathbf{S}_{\text{block}}$，然后执行子空间迭代：

1. 求解 block 系统得到新的子空间向量
2. M-正交化（双遍 Modified Gram-Schmidt）
3. Rayleigh-Ritz 投影求小特征值问题
4. 检查残差收敛

#### 方法二：SLEPc Krylov-Schur（`"slepc"`，需编译时启用）

将完整的 block 广义特征值问题：

$$\begin{bmatrix} \mathbf{A} & \mathbf{C} \\ \mathbf{0} & \mathbf{D} \end{bmatrix} \begin{bmatrix} \mathbf{e}_t \\ \tilde{\mathbf{e}}_z \end{bmatrix} = \beta^2 \begin{bmatrix} \mathbf{M} & \mathbf{0} \\ \mathbf{C}^T & \mathbf{0} \end{bmatrix} \begin{bmatrix} \mathbf{e}_t \\ \tilde{\mathbf{e}}_z \end{bmatrix}$$

转换为 PETSc 矩阵后交由 SLEPc 的 EPS 求解器处理，使用 Krylov-Schur 算法配合 shift-and-invert 谱变换。SLEPc 通常收敛更快（典型情况 ~14 次迭代 vs 内置求解器 ~160 次）。

SLEPc 求解器支持运行时参数覆盖，例如：
```bash
mpirun -np 1 ./build_slepc/eigensolve config.json -eps_monitor -eps_view
```

#### 方法三：ARPACK 隐式重启 Arnoldi（`"arpack"`，需编译时启用）

使用 PARPACK（并行 ARPACK）的隐式重启 Arnoldi 方法（Implicitly Restarted Arnoldi Method）求解归约后的特征值问题：

$$\mathbf{A}\, \mathbf{e}_t = \beta^2\, \mathbf{B}_{\text{eff}}\, \mathbf{e}_t$$

其中 $\mathbf{B}_{\text{eff}} = \mathbf{M} - \mathbf{C}\mathbf{D}^{-1}\mathbf{C}^T$。

采用 shift-and-invert 策略，定义算子 $\text{OP} = (\mathbf{A} - \sigma\mathbf{B}_{\text{eff}})^{-1}\mathbf{B}_{\text{eff}}$，通过 block 系统分解高效实现。ARPACK 通过反向通信（reverse communication）接口驱动迭代，每步需一次 SuperLU 直接求解和一次 D 矩阵求解。

### 边界条件

程序支持逐边独立设置边界条件（ymin、xmax、ymax、xmin 四个边界）：

| 类型 | 关键字 | 物理含义 | ND (E_t) | H1 (E_z) |
|------|--------|----------|----------|----------|
| **PEC** | `"pec"` | 理想电导体壁 (n×E=0) | 本质 BC (E_t=0) | 本质 BC (E_z=0) |
| **PMC** | `"pmc"` | 理想磁导体壁 (n×H=0) | 自然 BC | 自然 BC |
| **对称** | `"symmetric"` | 电壁对称面，等价于 PEC | 本质 BC | 本质 BC |
| **反对称** | `"antisymmetric"` | 磁壁对称面，等价于 PMC | 自然 BC | 自然 BC |
| **PML** | `"pml"` | 完美匹配层（吸收边界） | PEC (外层) | PEC (外层) |

#### PML（完美匹配层）

PML 通过实数坐标拉伸实现吸收，在 PML 区域引入拉伸函数：

$$s_x(x) = 1 + \sigma_{\max} \left(\frac{d}{d_{\text{pml}}}\right)^p$$

其中 $d$ 是进入 PML 的深度，$d_{\text{pml}}$ 是 PML 厚度，$p=2$，$\sigma_{\max}=5$。

PML 对各个双线性形式的修改：

| 双线性形式 | 标准形式 | PML 修改后 |
|-----------|---------|-----------|
| Curl-curl (ND) | $\int (\nabla_t \times \mathbf{e}_t)^2$ | $\int \frac{1}{s_x s_y} (\nabla_t \times \mathbf{e}_t)^2$ |
| Vector mass (ND) | $\int \varepsilon_r \mathbf{e}_t \cdot \mathbf{f}_t$ | $\int \varepsilon_r \Lambda_t \mathbf{e}_t \cdot \mathbf{f}_t$，$\Lambda_t = \text{diag}(s_y/s_x, s_x/s_y)$ |
| Diffusion (H1) | $\int \nabla e_z \cdot \nabla f_z$ | $\int \Lambda_t^{-1} \nabla e_z \cdot \nabla f_z$ |
| Scalar mass (H1) | $\int \varepsilon_r e_z f_z$ | $\int \varepsilon_r s_x s_y e_z f_z$ |

网格在 PML 方向自动延伸，PML 区域标记为 attribute 3，外层边界施加 PEC。

### 导模判据

有效折射率 $n_{\text{eff}} = \beta / k_0$ 满足 $n_{\text{clad}} < n_{\text{eff}} < n_{\text{core}}$ 的模式为导模。

## 依赖项

| 库 | 版本 | 用途 |
|---|------|------|
| MFEM | 4.9+ | 有限元框架（网格、空间、积分器、并行矩阵） |
| HYPRE | 3.x | 并行稀疏矩阵 (HypreParMatrix) |
| SuperLU_DIST | 9.x | Block 系统直接分解 |
| METIS/ParMETIS | 5.x | 网格分区、矩阵重排序 |
| LAPACK/OpenBLAS | — | 小规模 Rayleigh-Ritz 特征值问题 |
| nlohmann/json | 3.11+ | JSON 配置解析（CMake 自动下载） |
| MPI | — | 并行通信 |
| **SLEPc** | 3.24+ | **可选**：Krylov-Schur 特征值求解器 |
| **PETSc** | 3.24+ | **可选**：SLEPc 的依赖（并行线性代数） |
| **ARPACK-ng** | 3.9+ | **可选**：隐式重启 Arnoldi 特征值求解器 |

上述库通过 [Spack](https://spack.io/) 安装 MFEM 时自动包含。SLEPc/PETSc 和 ARPACK 为可选依赖，需单独安装：
```bash
spack install slepc       # SLEPc + PETSc
spack install arpack-ng   # ARPACK-ng (含 PARPACK)
```

## 构建

### 基础构建（内置求解器）

```bash
cd eigensolve
mkdir build && cd build
cmake ..
make -j4
```

CMake 会自动在 `~/spack/opt/spack/` 下搜索 MFEM 安装路径。也可手动指定：

```bash
cmake -DMFEM_DIR=/path/to/mfem/install ..
```

### 启用 ARPACK 求解器

```bash
mkdir build_arpack && cd build_arpack
cmake .. -DEIGENSOLVE_USE_ARPACK=ON
make -j4
```

CMake 自动检测 Spack 安装的 ARPACK-ng。也可手动指定路径：

```bash
cmake .. -DEIGENSOLVE_USE_ARPACK=ON -DARPACK_DIR=/path/to/arpack-ng
```

### 启用 SLEPc 求解器

```bash
mkdir build_slepc && cd build_slepc
cmake .. -DEIGENSOLVE_USE_SLEPC=ON
make -j4
```

CMake 自动检测 Spack 安装的 PETSc/SLEPc。也可手动指定路径：

```bash
cmake .. -DEIGENSOLVE_USE_SLEPC=ON \
         -DPETSC_DIR=/path/to/petsc \
         -DSLEPC_DIR=/path/to/slepc
```

## 使用方法

```bash
# 使用内置求解器（默认 PEC 边界）
mpirun -np 1 ./build/eigensolve eigensolve_config.json

# 不同边界条件
mpirun -np 1 ./build/eigensolve eigensolve_config_pec.json      # PEC
mpirun -np 1 ./build/eigensolve eigensolve_config_pml.json      # PML 吸收边界
mpirun -np 1 ./build/eigensolve eigensolve_config_sym.json      # 对称边界
mpirun -np 1 ./build/eigensolve eigensolve_config_antisym.json  # 反对称边界

# 使用 ARPACK 求解器（需要启用 ARPACK 构建）
mpirun -np 1 ./build_arpack/eigensolve eigensolve_config_arpack.json

# 使用 SLEPc 求解器（需要启用 SLEPc 构建）
mpirun -np 1 ./build_slepc/eigensolve eigensolve_config_slepc.json

# 多进程并行
mpirun -np 4 ./build/eigensolve eigensolve_config.json

# SLEPc 运行时参数（可覆盖默认设置）
mpirun -np 1 ./build_slepc/eigensolve config.json -eps_monitor -eps_view
```

## 配置文件说明

示例文件 `eigensolve_config.json`：

```json
{
    "wavelength": 1.55,
    "domain": {
        "width": 6.0,
        "height": 6.0
    },
    "waveguide": {
        "type": "rectangular",
        "width": 0.5,
        "height": 0.22,
        "radius": 0.5
    },
    "materials": {
        "core_index": 3.48,
        "cladding_index": 1.44
    },
    "boundary_conditions": {
        "ymin": "pec",
        "xmax": "pec",
        "ymax": "pec",
        "xmin": "pec",
        "pml_thickness": 1.0
    },
    "solver": {
        "order": 2,
        "num_modes": 6,
        "target_neff": 2.5,
        "refinements": 3,
        "element_type": "quad",
        "tolerance": 1e-4,
        "max_iterations": 200,
        "eigensolver": "builtin"
    },
    "output": {
        "directory": "results",
        "save_fields": true,
        "save_png": true
    }
}
```

### 参数详解

| 参数 | 类型 | 说明 |
|------|------|------|
| `wavelength` | float | 入射光波长（单位：μm） |
| `domain.width` | float | 计算域宽度（μm），应远大于波导截面 |
| `domain.height` | float | 计算域高度（μm） |
| `waveguide.type` | string | `"rectangular"` 或 `"circular"` |
| `waveguide.width` | float | 矩形波导宽度（μm） |
| `waveguide.height` | float | 矩形波导高度（μm） |
| `waveguide.radius` | float | 圆形波导半径（μm） |
| `materials.core_index` | float | 芯层折射率 |
| `materials.cladding_index` | float | 包层折射率 |
| `boundary_conditions.ymin` | string | 下边界条件：`"pec"`, `"pmc"`, `"pml"`, `"symmetric"`, `"antisymmetric"`（默认 `"pec"`） |
| `boundary_conditions.xmax` | string | 右边界条件（同上） |
| `boundary_conditions.ymax` | string | 上边界条件（同上） |
| `boundary_conditions.xmin` | string | 左边界条件（同上） |
| `boundary_conditions.pml_thickness` | float | PML 层厚度（μm），当任一边为 `"pml"` 时生效（默认 1.0） |
| `solver.order` | int | 有限元多项式阶数（建议 2-3） |
| `solver.num_modes` | int | 求解模式数 |
| `solver.target_neff` | float | 目标有效折射率（shift 参数，取芯层和包层折射率之间的值） |
| `solver.refinements` | int | 芯层边界自适应加密次数 |
| `solver.element_type` | string | `"triangle"` 或 `"quad"` |
| `solver.tolerance` | float | 特征值残差收敛容限 |
| `solver.max_iterations` | int | 最大迭代次数 |
| `solver.eigensolver` | string | `"builtin"`（默认）、`"slepc"` 或 `"arpack"`，选择特征值求解器 |
| `output.directory` | string | 场输出目录 |
| `output.save_fields` | bool | 是否保存场分布 VTU 文件 |
| `output.save_png` | bool | 是否生成 PNG 场分布图 |

## 输出

### 控制台输出

程序打印每个导模的传播常数和有效折射率：

```
╔═══════════════════════════════════════════════════════╗
║           Waveguide Mode Analysis Results             ║
╠═══════╦═══════════════╦═══════════════╦═══════════════╣
║ Mode  ║     n_eff     ║    beta       ║   beta^2      ║
║       ║               ║  (1/um)       ║  (1/um^2)     ║
╠═══════╬═══════════════╬═══════════════╬═══════════════╣
║     1 ║    2.44650183 ║      9.917306 ║       98.3530 ║
║     2 ║    1.79494808 ║      7.276123 ║       52.9420 ║
║     3 ║    1.48252560 ║      6.009666 ║       36.1161 ║
╚═══════╩═══════════════╩═══════════════╩═══════════════╝
```

### 场分布文件

每个模式保存为 ParaView VTU 格式和 PNG 图片：

```
results/
├── mode_0/              # TE₀ 基模（VTU 格式）
│   ├── mode_0.pvd
│   └── Cycle000000/
│       ├── data.pvtu
│       └── proc000000.vtu
├── mode_0.dat           # 场数据（均匀网格采样）
├── mode_0_Ex.png        # Ex 分量分布图
├── mode_0_Ey.png        # Ey 分量分布图
├── mode_0_abs.png       # |Et| 模值分布图
├── mode_1/
│   └── ...
└── mode_N/
    └── ...
```

- 使用 [ParaView](https://www.paraview.org/) 打开 `.pvd` 文件即可可视化横向电场 $\mathbf{E}_t$ 分布
- PNG 图片由 `plot_modes.py` 自动生成（需安装 numpy、matplotlib），也可手动运行：
  ```bash
  python3 plot_modes.py results
  ```

## 项目结构

```
eigensolve/
├── CMakeLists.txt                  # 构建系统（支持 -DEIGENSOLVE_USE_SLEPC=ON）
├── eigensolve_config.json          # 示例配置（内置求解器，默认 PEC）
├── eigensolve_config_pec.json     # PEC 边界配置
├── eigensolve_config_pml.json     # PML 吸收边界配置
├── eigensolve_config_sym.json     # 对称边界配置
├── eigensolve_config_antisym.json # 反对称边界配置
├── eigensolve_config_arpack.json  # ARPACK 求解器配置
├── eigensolve_config_slepc.json   # SLEPc 求解器配置
├── plot_modes.py                   # PNG 场分布图生成脚本
├── src/
│   ├── main.cpp                    # 入口：MPI/SLEPc 初始化、配置解析、主流程
│   ├── mesh_generator.hpp/cpp      # 2D 网格生成 + PML 扩展 + 自适应加密 + 材料标记
│   ├── waveguide_solver.hpp/cpp    # FE 空间、PML 系数类、矩阵组装、边界条件处理
│   ├── eigensolver.hpp/cpp         # EigensolverBase 接口 + 内置子空间迭代
│   ├── slepc_eigensolver.hpp/cpp   # SLEPc 特征值求解器（可选编译）
│   └── arpack_eigensolver.hpp/cpp # ARPACK 特征值求解器（可选编译）
├── config.json                     # （旧 Palace 配置，供参考）
└── waveguide.geo                   # （旧 Gmsh 几何，供参考）
```

## 验证基准

标准硅基条形波导（500 nm × 220 nm，SiO₂ 包层，λ = 1550 nm）：

### 不同边界条件对比

| 模式 | PEC n_eff | PML n_eff | Symmetric n_eff | Antisymmetric n_eff |
|------|-----------|-----------|-----------------|---------------------|
| 1 | 2.44646 | 2.50514 | 2.44646 | 2.44646 |
| 2 | 1.79489 | 2.50514 | 1.79489 | 1.79489 |
| 3 | 1.48245 | 2.50514 | 1.48245 | 1.48243 |

- **PML 边界**的模式 n_eff 值更高（~2.505），因为 PML 吸收了辐射模，改善了导模精度
- **PEC** 与**对称/反对称**边界在全域计算中结果相近（对称=电壁，反对称=磁壁）
- 对称/反对称边界可用于利用波导的几何对称性，减小计算域为 1/4 或 1/2

### 求解器对比

| 模式 | 内置求解器 n_eff | SLEPc n_eff | ARPACK n_eff | 文献参考值 |
|------|---------------|------------|-------------|-----------|
| TE₀ | 2.44646 | 2.44646 | 2.44646 | ~2.4 |
| TM₀ | 1.79489 | 1.79489 | 1.79489 | ~1.7–1.8 |

三种求解器结果完全一致（差异 < 1e-5）。

| 求解器 | Arnoldi/子空间迭代次数 | OP*x 调用次数 | 备注 |
|--------|---------------------|-------------|------|
| 内置 (builtin) | ~160 次 | ~160 次 | 简单子空间迭代 |
| SLEPc | ~14 次 | — | Krylov-Schur |
| ARPACK | ~14 次 | ~163 次 | 隐式重启 Arnoldi |

## 使用建议

- **计算域大小**：建议为波导截面尺寸的 10 倍以上，PEC 边界需更大域以减少边界反射影响
- **边界条件选择**：
  - **PEC**（默认）：最简单，但域边界处存在反射，需较大计算域
  - **PML**：推荐用于开放波导，通过吸收层消除边界反射，结果更准确
  - **对称/反对称**：利用波导几何对称性，可将计算域缩小为 1/2 或 1/4，大幅降低计算量
- **PML 参数**：`pml_thickness` 建议设为 0.5–2.0 μm，过薄吸收不足，过厚增加计算量
- **target_neff**：设为预期有效折射率附近的值，通常取 $(n_{\text{core}} + n_{\text{clad}}) / 2$ 即可
- **refinements**：3 次加密通常足够，增加到 4-5 可提高精度但增加计算量
- **order**：2 阶通常性价比最优；3 阶可用于高精度需求
- **element_type**：三角形适合一般情况；四边形在规则区域效率更高

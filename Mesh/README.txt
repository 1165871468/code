批量测线网格生成

1. 将 build_topo_mesh.py 与本包 read_file.py 放在同一目录。
2. 修改主脚本顶部 CONFIG_FILE、LINE_FILE、OUTPUT_DIR 三个路径。
   默认相对路径基于运行时工作目录。
3. 运行 python build_topo_mesh.py，需要 NumPy。

配置读取调用 read_file.read_forward_config：
  dis：物理区水平单元间隔，单位 m。
  Layx：左右各增加的 PML 单元数。
  Layz：底部增加的 PML 单元数，无顶部 PML。
  model_size：参考长度、深度两个数，只使用第二个数作为物理区深度。
  原配置的末行英文说明含 horizontal and vertical scales 时，也可识别。
  推荐明确写成：100000 40000 c: Model size // model_size

line.txt 必须是原来的九列格式：
source_id receiver_id xs ys zs xe ye ze length
第九列长度作为模型长度依据，不使用源台距离，也不修改测线方向。
模型物理区长度 = ceil(length/dis)*dis；x 从 0 起，扩展部分加在右端。
深度保持配置值；nz=ceil(depth/dis)，基准 dz=depth/nz。
PML 再向左右及底部扩展，不包含在物理区尺寸内。
正弦地形参数 AMPLITUDE、HALF_WAVES 保留在脚本顶部，默认振幅为 0。

输出：
Mesh_s1_r1.txt 等：保留原四节点网格格式、节点编号和材料号。
材料号 1 为物理区，2 为 PML。
mesh_summary.txt：带 # 表头，记录每条测线的尺寸和网格统计。
WRITE_INTERMEDIATE=True 可输出带源台编号的三个中间文件。
同名输出文件在再次运行时会覆盖；重复源台编号会报错。

read_file.py 只新增了无 // 标记的模型尺寸行识别分支，其余接口保持不变。
本程序以此补充版为配套文件；不要使用不能识别深度字段的旧版读取模块。

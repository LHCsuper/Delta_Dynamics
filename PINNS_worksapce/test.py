import sys

print("==== Python 环境信息 ====")
print("Python 可执行文件路径:", sys.executable)
print("Python 版本:", sys.version)

print("\n==== PyTorch 检测 ====")
try:
    import torch
    print("PyTorch 已安装")
    print("PyTorch 版本:", torch.__version__)
    print("CUDA 是否可用:", torch.cuda.is_available())

    if torch.cuda.is_available():
        print("CUDA 版本:", torch.version.cuda)
        print("GPU 数量:", torch.cuda.device_count())
        print("当前 GPU 名称:", torch.cuda.get_device_name(0))
    else:
        print("当前为 CPU 环境，或 CUDA 未配置成功")

    x = torch.tensor([1.0, 2.0, 3.0])
    y = x * 2
    print("张量测试结果:", y)

except ImportError as e:
    print("未检测到 PyTorch")
    print("报错信息:", e)
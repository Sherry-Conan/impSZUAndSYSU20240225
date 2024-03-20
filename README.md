# 401中山大学深圳大学深圳科技大学联合实验
于2024年2月25日上传
## 需要编译后运行
## 编译环境要求
### root6.30或以上
### C++17或以上
### GCC8或以上
### Clang6或以上

## 作者建议
推荐使用archlinux 虚拟机和双系统教程 https://arch.icekylin.online/ wsl教程(不建议安装桌面系统KDE) https://gitee.com/regentsai/wsl_arch_kde#%E7%BB%88%E7%AB%AF%E6%93%8D%E4%BD%9Cwsl%E7%9A%84%E5%9F%BA%E6%9C%AC%E7%94%A8%E6%B3%95

于2024年3月20日
### 提交了比较大的改动
## 多出的ipynb文件不一定好用，需要自己改，主要是我手里没有本次实验的数据，只好用之前的了。
## 对其他文件也进行了改动，我个人喜欢在后处理中引入时间对齐，若想改动需要：
### 在 imp.h 文件中的73行变 TimeLa[LaNum - 1] 为 TimeLa[LaNum]
### 在 imp.cpp 文件中的74行 添加 + TimeLa[Ch[i]]

# 以后如有改动会再更新
## 推荐使用archlinux 虚拟机和双系统教程 https://arch.icekylin.online/ wsl教程(不建议安装桌面系统KDE) https://gitee.com/regentsai/wsl_arch_kde#%E7%BB%88%E7%AB%AF%E6%93%8D%E4%BD%9Cwsl%E7%9A%84%E5%9F%BA%E6%9C%AC%E7%94%A8%E6%B3%95

# コンパイラとフラグの設定
CXX = g++
CXXFLAGS = -O3 -w -std=c++17
OPENCV_FLAGS = `pkg-config --cflags --libs opencv4`
GLUT_FLAGS = -I/usr/local/Cellar/freeglut/3.6.0/include -L/usr/local/Cellar/freeglut/3.6.0/lib -lglut -framework OpenGL
EIGEN_FLAGS = -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3
BOOST_FLAGS = -I/usr/local/opt/boost/include

# OpenMPの設定
ifeq ($(shell uname), Darwin) # macOSの場合
    OPENMP_FLAGS = -Xpreprocessor -fopenmp -I/usr/local/Cellar/libomp/18.1.8/include
    OPENMP_LIBS = -L/usr/local/Cellar/libomp/18.1.8/lib -lomp
else # Linuxの場合
    OPENMP_FLAGS = -fopenmp
    OPENMP_LIBS = -lgomp
endif

# 出力ファイルとオブジェクトファイル
TARGET = Bin/membrane
OBJS = Obj/main.o Obj/Grid.o Obj/Particle_Cloud.o Obj/Particle.o Obj/Grid_Point.o

# メインターゲット
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OPENCV_FLAGS) $(GLUT_FLAGS) $(OPENMP_LIBS) $(BOOST_FLAGS)

# オブジェクトファイルの生成ルール
Obj/main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(OPENCV_FLAGS) $(EIGEN_FLAGS) $(GLUT_FLAGS) $(OPENMP_FLAGS) $(BOOST_FLAGS)

Obj/Grid.o: Source/Grid.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(OPENCV_FLAGS) $(EIGEN_FLAGS) $(GLUT_FLAGS) $(OPENMP_FLAGS) $(BOOST_FLAGS)

Obj/Grid_Point.o: Source/Grid_Point.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(GLUT_FLAGS) $(OPENMP_FLAGS) $(BOOST_FLAGS)

Obj/Particle_Cloud.o: Source/Particle_Cloud.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(GLUT_FLAGS) $(OPENMP_FLAGS) $(BOOST_FLAGS)

Obj/Particle.o: Source/Particle.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(GLUT_FLAGS) $(OPENMP_FLAGS) $(BOOST_FLAGS)

# クリーンアップ
clean:
	rm -f $(OBJS) $(TARGET)

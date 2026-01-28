
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <windows.h>
#include <random>
#include <format>

using namespace std;


class TreeNode
{
public:
	vector<int> v;
	vector<vector<int>> genSystem;
	vector<TreeNode> Mv;
	int level;
	int maxLevel;
};

class TensorTreeNode
{
public:
	TreeNode v1;
	TreeNode v2;
	vector<TensorTreeNode> Mv;
	int value;
	int level;
	int maxLevel;
};



class ReedMallerTreeBuilder{

public:
	// декодирующее дерево для кода Рида-Маллера
	TreeNode tree;
	// параметры кода Рида-Маллера
	int r, m;

	// конструктор класса ReedMallerTreeBuilder
	ReedMallerTreeBuilder(int R, int M)
	{
		r = R;
		m = M;
		allElements = createMap(m);
		tree = makeRMWBTree(r, m);
	}

	// вывод декодирующего дерева
	void printTree(const string& name)
	{
		cout << name << ":" << endl << endl;
		printTreeRec(tree);
		cout << "--------------------------------------" << endl;
	}

private:

	// структура, задающая порядок на группе
	map<vector<int>, int> allElements;

	// покоординатное сложение двух векторов
	vector<int> sumVectors(vector<int> numbers, vector<vector<int>> vectors)
	{
		vector<int> res(vectors[0].size(), 0);

		for (int i = 0; i < res.size(); ++i)
		{
			int s = 0;
			for (int j = 0; j < vectors.size(); ++j)
			if (numbers[j])
			{
				s += vectors[j][i];
			}
			res[i] = s % 2;
		}
		return res;
	}

	// вычисление списка всех векторов, порожденных векторами genSystem
	vector<vector<int>> allWords(vector<vector<int>> genSystem)
	{
		vector<vector<int>> words;
		vector<int> v(genSystem.size(), 0);
		bool flag = true;
		while (flag)
		{
			int k = -1;
			for (int i = 0; i < genSystem.size(); ++i)
			if (v[i] == 0)
			{
				k = i;
				v[i] = 1;
				for (int j = k - 1; j >= 0; --j)
					v[j] = 0;

				vector<int> vec = sumVectors(v, genSystem);
				words.push_back(vec);
				break;
			}
			if (k == -1)
				flag = false;
		}
		return words;
	}

	// алгоритм Гаусса приведения матрицы к ступенчатому виду
	void gauss(vector<vector<int>> &matr, vector<int> &errorCols)
	{
		for (int i = 0, j = 0; i < matr.size() && j < matr[0].size(); ++i, ++j)
		{
			if (matr[i][j] == 0)
			{
				int num = -1;
				for (int k = i; k < matr.size(); ++k)
				{
					if (matr[k][j] != 0)
					{
						num = k;
						break;
					}
				}
				if (num != -1)
					swap(matr[i], matr[num]);
				else
				{
					errorCols[j] = 1;
					--i;
					continue;
				}
			}
			else
			{
				for (int k = i + 1; k < matr.size(); ++k)
				{
					if (matr[k][i] == 1)
					for (int l = i; l < matr[i].size(); ++l)
						matr[k][l] = (matr[k][l] + matr[i][l]) % 2;
				}
			}
		}
	}

	// добавление одного уровня декодирующего дерева для кода Рида-Маллера
	void addLevel(vector<TreeNode*> &V)
	{
		vector<TreeNode*> newV;
		int size = pow(2, m);
		for (int i = 0; i < V.size(); ++i)
		{
			vector<vector<int>> oldGenSystem = V[i]->genSystem;
			// система порождающих элемментов для нового узла дерева
			vector<vector<int>> newGenSystem;
			vector<int> errs(m, 0);

			// алгоритм Гаусса приведения матрицы к ступенчатому виду
			gauss(oldGenSystem, errs);

			for (int j = 0; j < errs.size(); ++j)
			if (errs[j] == 1)
			{
				vector<int> v(errs.size(), 0);
				v[j] = 1;
				newGenSystem.push_back(v);
			}

			int num = -1;
			for (int j = 0; j < oldGenSystem[oldGenSystem.size() - 1].size(); ++j)
			if (oldGenSystem[oldGenSystem.size() - 1][j] == 1)
			{
				num = j;
				break;
			}


			if (num != oldGenSystem[oldGenSystem.size() - 1].size() - 1)
			for (int j = num + 1; j < oldGenSystem[oldGenSystem.size() - 1].size(); ++j)
			{
				vector<int> v(errs.size(), 0);
				v[j] = 1;
				newGenSystem.push_back(v);
			}

			// вычисление списка всех векторов, порожденных векторами newGenSystem
			newGenSystem = allWords(newGenSystem);

			// построение M-ортогональной системы для текущего узла дерева
			for (int j = 0; j < newGenSystem.size(); ++j)
			{
				TreeNode ti;
				ti.genSystem = V[i]->genSystem;
				ti.genSystem.push_back(newGenSystem[j]);
				ti.level = V[i]->level + 1;
				ti.maxLevel = V[i]->maxLevel;

				vector<int> vi(size, 0);
				vi[0] = 1;
				vector<vector<int>> words = allWords(ti.genSystem);
				for (int i = 0; i < words.size(); ++i)
					vi[allElements[words[i]]] = 1;

				ti.v = vi;

				V[i]->Mv.push_back(ti);
			}

			for (int j = 0; j < V[i]->Mv.size(); ++j)
			{
				TreeNode* Vi = &V[i]->Mv[j];
				newV.push_back(Vi);
			}
		}

		V = newV;
	}

	// построение декодирующего дерева для RM(R,M)-кода Рида-Маллера
	TreeNode makeRMWBTree(int R, int M)
	{
		int N = pow(2, M);
		TreeNode tree;
		vector<int> v(N, 0);
		v[0] = 1;
		tree.v = v;
		tree.level = 0;
		tree.maxLevel = R + 1;

		// построение корня декодирующего дерева для кода Рида-Маллера
		for (auto it = allElements.begin(); it != allElements.end(); ++it)
		{
			TreeNode ti;
			vector<int> vi(N, 0);

			ti.genSystem.push_back(it->first);

			vi[0] = 1;
			vector<vector<int>> words = allWords(ti.genSystem);
			for (int i = 0; i < words.size(); ++i)
				vi[allElements[words[i]]] = 1;

			ti.v = vi;
			ti.level = 1;
			ti.maxLevel = tree.maxLevel;
			tree.Mv.push_back(ti);
		}

		// список всех узлов дерева, находящихся на одном уровне
		vector<TreeNode*> V;
		for (int i = 0; i < tree.Mv.size(); ++i)
		{
			TreeNode* Vi = &tree.Mv[i];
			V.push_back(Vi);
		}

		for (int i = 1; i < tree.maxLevel; ++i)
			// добавление одного уровня декодирующего дерева
			addLevel(V);
		return tree;
	}

	// функция пстроения вспомогательной структуры allElements
	map<vector<int>, int> createMap(int n)
	{
		map<vector<int>, int> allElements;
		vector<int> v(n, 0);
		bool flag = true;
		int l = 1;
		while (flag)
		{
			int k = -1;
			for (int i = 0; i < n; ++i)
			if (v[i] == 0)
			{
				k = i;
				v[i] = 1;
				for (int j = k - 1; j >= 0; --j)
					v[j] = 0;

				allElements.insert(pair<vector<int>, int>(v, l));
				l++;
				break;
			}
			if (k == -1)
				flag = false;
		}
		return allElements;
	}

	// вывод декодирующего дерева
	void printTreeRec(TreeNode t)
	{
		vector<int> val = t.v;
		for (int j = 0; j < t.level * 2; ++j)
			cout << "  ";
		for (int i = 0; i < val.size(); ++i)
		{
			cout << val[i];
		}
		cout << endl;

		if (t.Mv.size() != 0)
		for (int i = 0; i < t.Mv.size(); ++i)
			printTreeRec(t.Mv[i]);
	}
};



class TensorDecoder{

public:
	// конструктор класса TensorDecoder
	// на вход подаются декодирующие деревья для производных кодов
	TensorDecoder(TreeNode a, TreeNode b)
	{
		TensorTreeNode tensorTree = buildTensorTree(a, b);
		tensorTrees = cloneTree(tensorTree);
	}

	// функция, осуществляющая декодирование вектора x
	vector<int> decode(vector<int> x)
	{
		for (int i = 0; i < x.size(); ++i)
		{
			x[i] += decodeTensorBit(tensorTrees[i], x);
			x[i] %= 2;
		}
		return x;
	}

	// вывод декодирующего дерева
	void printTree(int n, const string& name)
	{
		cout << name << ":" << endl << endl;
		printTensorTreeRec(tensorTrees[n]);
		cout << "--------------------------------------" << endl;
	}

private:

	// список декодирующих деревьев
	vector<TensorTreeNode> tensorTrees;

	// функция, возвращающая тензорное произведение векторов a и b
	static vector<int> tensorProduct(vector<int> a, vector<int> b)
	{
		int k1 = 0, k2 = 0;
		vector<int> v(a.size()*b.size(), 0);
		for (int i = 0; i < a.size(); ++i)
		for (int j = 0; j < b.size(); ++j){
			v[i*b.size() + j] = a[i] * b[j];
		}
		return v;
	}

	// построение M-ортогональной системы для тензороного произведения векторов a и b
	vector<TensorTreeNode> mOrth(TreeNode a, TreeNode b, int level, int maxLevel)
	{
		vector<TensorTreeNode> Mv;
		for (int i = 0; i < b.Mv.size(); ++i)
		{
			TensorTreeNode ti;
			ti.v1 = a;
			ti.v2 = b.Mv[i];
			ti.level = level;
			ti.maxLevel = maxLevel;
			Mv.push_back(ti);
		}
		for (int i = 0; i < a.Mv.size(); ++i)
		{
			TensorTreeNode ti;
			ti.v1 = a.Mv[i];
			ti.v2 = b;
			ti.level = level;
			ti.maxLevel = maxLevel;
			Mv.push_back(ti);
		}
		return Mv;
	}

	// функция построения вспомогательного декодирующего дерева
	TensorTreeNode buildTensorTree(TreeNode a, TreeNode b)
	{
		TensorTreeNode t;
		t.v1 = a; t.v2 = b;
		t.level = 0;
		t.maxLevel = a.maxLevel + b.maxLevel - 1;
		t.Mv = mOrth(a, b, 1, t.maxLevel);
		// массив, который хранит список всех узлов 
		// одного уровня декодирующего дерева
		vector<TensorTreeNode*> V;
		for (int i = 0; i < t.Mv.size(); ++i)
		{
			TensorTreeNode *vi = &t.Mv[i];
			V.push_back(vi);
		}

		for (int i = 1; i < a.maxLevel + b.maxLevel - 1; ++i)
		{
			vector<TensorTreeNode*> Vnew;
			for (int j = 0; j < V.size(); ++j)
			if (V[j]->v1.level < V[j]->v1.maxLevel  && V[j]->v2.level < V[j]->v2.maxLevel)
			{
				V[j]->Mv = mOrth(V[j]->v1, V[j]->v2, i + 1, t.maxLevel);
				for (int k = 0; k < V[j]->Mv.size(); ++k)
				{
					TensorTreeNode *vi = &(V[j]->Mv[k]);
					Vnew.push_back(vi);
				}
			}
			V = Vnew;
		}

		return t;
	}

	// функция, вычисляющая скалярное произведение векторов a и b
	int innerProduct(vector<int> a, vector<int> b, int mod)
	{
		int s = 0;
		for (int i = 0; i < a.size(); ++i)
			s += a[i] * b[i];
		return s % mod;
	}

	// функция, вычисляющая значения дополнительных меток, необходимых для декодирования
	vector<int> addVals(TensorTreeNode tree, vector<int> x)
	{
		vector<int> l;
		if (tree.level == tree.v1.maxLevel + tree.v2.maxLevel - 2)
		{
			for (int i = 0; i < tree.v2.Mv.size(); ++i)
			for (int j = tree.v2.Mv.size(); j < tree.v1.Mv.size() + tree.v2.Mv.size(); ++j)
				l.push_back((tree.Mv[i].value + tree.Mv[j].value + innerProduct(tensorProduct(tree.Mv[j].v1.v, tree.Mv[i].v2.v), x, 2)) % 2);
		}
		else
		{
			for (int i = 0; i < tree.v2.Mv.size(); ++i)
			for (int j = 0; j < tree.v1.Mv.size(); ++j)
			{
				// поиск значения метки для третьего слагаемого конструкции 3 леммы 1
				int k = 0;
				if (tree.Mv[tree.v2.Mv.size() + j].Mv.size() == 0)
					k = tree.Mv[i].Mv[tree.Mv[i].v2.Mv.size() + j].value;
				else
					k = tree.Mv[tree.v2.Mv.size() + j].Mv[i].value;
				l.push_back((tree.Mv[i].value + tree.Mv[tree.v2.Mv.size() + j].value + k) % 2);
			}
		}

		return l;
	}

	// функция, осуществляющая процедуру мажоритарного голосования
	int majorVote(vector<int> l, int mod)
	{
		l.push_back(0);
		vector<int> v(mod);
		for (int i = 0; i < l.size(); ++i)
			v[l[i]]++;

		int max = -1, k = -1;
		for (int i = 0; i < v.size(); ++i)
		if (v[i] > max)
		{
			max = v[i];
			k = i;
		}

		return k;
	}

	// функция, осуществляющая декодирование координаты входного вектора x
	int decodeTensorBit(TensorTreeNode &tree, vector<int> x)
	{
		if (tree.Mv.size() == 0)
		{
			// вычисление скалярного произведения
			// для узлов дерева, принадлежащих ортогональному коду
			tree.value = innerProduct(x, tensorProduct(tree.v1.v, tree.v2.v), 2);
		}
		else
		{
			vector<int> l;
			for (int i = 0; i < tree.Mv.size(); ++i)
			{
				l.push_back(decodeTensorBit(tree.Mv[i], x));
			}

			// вычисление дополнительных меток для декодирования текущего узла
			vector<int> additionalValues = addVals(tree, x);
			l.insert(l.end(), additionalValues.begin(), additionalValues.end());
			// вычисление метки для текущего узла путем мажоритарного голосования
			tree.value = majorVote(l, 2);
		}

		return tree.value;
	}

	// применение элемента группы с номером g к вектору v
	vector<int> applyGroupElement(vector<int> v, int g) {

		int n = v.size();
		int m = log2(n);
		vector<vector<int>> numbers;
		numbers.reserve(n);

		// получаем множество всех вектор-номеров (элементов группы)
		for (int val = 0; val < n; ++val) {
				vector<int> num(m + 1);
				// заполняем двоичное представление (старший бит слева)
				for (int i = 0; i < m; ++i) {
					int bitPos = m - 1 - i;       // позиция бита от старшего к младшему
					num[i] = (val >> bitPos) & 1;   // 0 или 1
				}
				num[m] = val; // последний элемент — десятичное значение
				numbers.push_back(num);
			}

		// получаем вектор-номер элемента g
		vector<int> ge = numbers[g];

		// применяем элемент g ко всем вектор-номерам
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				numbers[i][j] = (numbers[i][j] + ge[j]) % 2;

		// упорядочиваем по вектор-номерам
		sort(numbers.begin(), numbers.end(),
              [](const vector<int>& x, const vector<int>& y) {
                  int n = (int)x.size();
                  for (int i = 0; i < n - 1; ++i) {
                      if (x[i] < y[i]) return true;
                      if (x[i] > y[i]) return false;
                  }
                  return false;
              });

		vector<int> res(n);
		for (int i = 0; i < n; ++i)
			res[i] = v[numbers[i][numbers[i].size()-1]];

		return res;
	}

	// применения ко всем узлам дерева элемента s1 по первому аргументу и s2 по второму
	void shiftTree(TensorTreeNode &tree, int s1, int s2)
	{

		tree.v1.v = applyGroupElement(tree.v1.v, s1);
		tree.v2.v = applyGroupElement(tree.v2.v, s2);

		if (tree.Mv.size() == 0)
			return;

		for (int i = 0; i < tree.Mv.size(); ++i)
		{
			shiftTree(tree.Mv[i], s1, s2);
		}
	}

	// функция построения декодирующих деревьев для каждой координаты
	vector<TensorTreeNode> cloneTree(TensorTreeNode tree)
	{
		vector<TensorTreeNode> allTrees;
		for (int i = 0; i < tree.v1.v.size(); ++i)
		for (int j = 0; j < tree.v2.v.size(); ++j)
		{
			TensorTreeNode shifted = tree;
			// применение ко всем узлам дерева элемента i по первому аргументу и j по второму
			shiftTree(shifted, i, j);
			allTrees.push_back(shifted);
		}

		return allTrees;
	}

	// вывод декодирующего дерева
	void printTensorTreeRec(TensorTreeNode t)
	{
		vector<int> val = tensorProduct(t.v1.v, t.v2.v);
		for (int j = 0; j < t.level * 2; ++j)
			cout << "  ";
		for (int i = 0; i < val.size(); ++i)
		{
			cout << val[i];
		}
		cout << "  " << endl;

		if (t.Mv.size() != 0)
		for (int i = 0; i < t.Mv.size(); ++i)
			printTensorTreeRec(t.Mv[i]);

	}
};



class TensorEncoder {
public:
    // построение порождающей матрицы кода RM(r, m)
    // возвращает матрицу размера k x n
    vector<vector<int>> buildRMGeneratorMatrix(int r, int m) {

        int n = 1 << m;
        vector<vector<int>> points(n, vector<int>(m, 0));
        for (int idx = 0; idx < n; ++idx) {
            for (int j = 0; j < m; ++j) {
                points[idx][j] = (idx >> j) & 1;
            }
        }


        vector<vector<int>> monomials;
        int total = 1 << m;
        for (int mask = 0; mask < total; ++mask) {
            int wt = __builtin_popcount(mask);
            if (wt > r) continue;
            vector<int> a(m, 0);
            for (int j = 0; j < m; ++j) {
                a[j] = (mask >> j) & 1;
            }
            monomials.push_back(a);
        }

        int k = (int)monomials.size();
        vector<vector<int>> G(k, vector<int>(n, 0));


        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < n; ++j) {
                int val = 1;
                for (int t = 0; t < m; ++t) {
                    if (monomials[i][t] == 1) {
                        val &= points[j][t]; 
                    }
                }
                G[i][j] = val & 1;
            }
        }

        return G;
    }

    // вычисление матрицы тензорного произведения двух входных матриц
    vector<vector<int>> buildTensorMatrix(const vector<vector<int>>& matr1, const vector<vector<int>>& matr2) {
        int k1 = (int)matr1.size();
        int n1 = k1 ? (int)matr1[0].size() : 0;
        int k2 = (int)matr2.size();
        int n2 = k2 ? (int)matr2[0].size() : 0;

        vector<vector<int>> res(k1 * k2, vector<int>(n1 * n2, 0));

        for (int i = 0; i < k1; ++i) {
            for (int j = 0; j < k2; ++j) {
                int row = i * k2 + j;
                for (int a = 0; a < n1; ++a) {
                    for (int b = 0; b < n2; ++b) {
                        res[row][a * n2 + b] =
                            (matr1[i][a] * matr2[j][b]) & 1;
                    }
                }
            }
        }

        return res;
    }

    // умножение вектора v (длина k) на матрицу gen_matr (k x n)
    vector<int> encodeTensor(const vector<int>& v, const vector<vector<int>>& gen_matr) {
        int k = (int)gen_matr.size();
        if (k == 0) return {};
        int n = (int)gen_matr[0].size();

        vector<int> res(n, 0);
        for (int i = 0; i < k; ++i) {
            if (v[i] == 0) continue;
            for (int j = 0; j < n; ++j) {
                res[j] ^= gen_matr[i][j];
            }
        }
        return res;
    }

	// сумма(разница) двух векторов
	vector<int> vectorDiff(const vector<int>& a, const vector<int>& b) {
		vector<int> diff;

		if (a.size() != b.size()) {
			return diff;
		}

		diff.resize(a.size());
		for (size_t i = 0; i < a.size(); ++i) {
			diff[i] = (a[i] + b[i]) % 2;
		}

    	return diff;
	}

	// генерация случайного информационного вектора длины k
    vector<int> generateRandomInfoVector(int k) {
        vector<int> v(k, 0);
        static random_device rd;
        static mt19937 gen(rd());
        uniform_int_distribution<int> bitDist(0, 1);
        for (int i = 0; i < k; ++i) {
            v[i] = bitDist(gen);
        }
        return v;
    }

    // добавление к вектору v случайной ошибки веса t, возвращает новый вектор
    vector<int> addNoise(vector<int> v, int t) {
        int n = (int)v.size();
        if (t <= 0 || n == 0) return v;
        if (t > n) t = n;

        static random_device rd;
        static mt19937 gen(rd());

        vector<int> pos(n);
        for (int i = 0; i < n; ++i) pos[i] = i;
        shuffle(pos.begin(), pos.end(), gen);

        for (int i = 0; i < t; ++i) {
            int p = pos[i];
            v[p] ^= 1;
        }
        return v;
    }

	// вывод матрицы в консоль
    void printMatrix(const vector<vector<int>>& matr, const string& name) {
        cout << name << " (" << matr.size() << " x "
             << (matr.empty() ? 0 : (int)matr[0].size()) << "):" << endl << endl;
        for (int i = 0; i < (int)matr.size(); ++i) {
            for (int j = 0; j < (int)matr[i].size(); ++j) {
                cout << matr[i][j];
            }
            cout << endl;
        }
        cout << "--------------------------------------" << endl;
    }

    // вывод вектора в консоль
    void printVector(const vector<int>& v, const string& name) {
        cout << name << " (length = " << v.size() << "):" << endl << endl;
        for (int i = 0; i < (int)v.size(); ++i) {
            cout << v[i];
        }
        cout << endl << "--------------------------------------" << endl;
    }

};



// кодирование и декодирование случайного информационного вектора с помощью тензорного произведения двух кодов Рида-Маллера
void runDecoderRand(int r1, int m1, int r2, int m2, int t) {

	ReedMallerTreeBuilder tb1(r1, m1);		// декодирующее дерево для кода Рида-Маллера 1
	ReedMallerTreeBuilder tb2(r2, m2);		// декодирующее дерево для кода Рида-Маллера 2
	TensorDecoder td(tb1.tree, tb2.tree);	// декодирующее дерево для тензорного произведения двух кодов Рида-Маллера

	cout << format("r1={}, m1={}, r2={}, m2={}, t={}", r1, m1, r2, m2, t) << endl;

	tb1.printTree(format("Декодирующее дерево для первой координаты кода RM({},{})", r1, m1));
	tb2.printTree(format("Декодирующее дерево для первой координаты кода RM({},{})", r2, m2));
	td.printTree(0, format("Вспомогательное декодирующее дерево для первой координаты тензорного произведения RM({},{}) и RM({},{})", r1, m1, r2, m2));

	TensorEncoder enc;

	vector<vector<int>> G1 = enc.buildRMGeneratorMatrix(r1, m1);	// порождающая матрица для кода Рида-Маллера 1
	vector<vector<int>> G2 = enc.buildRMGeneratorMatrix(r2, m2);	// порождающая матрица для кода Рида-Маллера 2
	vector<vector<int>> Gt = enc.buildTensorMatrix(G1, G2);			// порождающая матрица для тензорного произведения двух кодов Рида-Маллера

	enc.printMatrix(G1, format("Порождающая матрица G1 кода RM({},{})", r1, m1));
    enc.printMatrix(G2, format("Порождающая матрица G2 кода RM({},{})", r2, m2));
    enc.printMatrix(Gt, format("Порождающая матрица Gt тензорного произведения RM({},{}) и RM({},{})", r1, m1, r2, m2));

	int k = (int)Gt.size(); // длина информационного вектора

	vector<int> u = enc.generateRandomInfoVector(k);	// случайный информационный вектор
	vector<int> c = enc.encodeTensor(u, Gt);     		// кодирование вектора u
	vector<int> y = enc.addNoise(c, t);          		// добавили ошибку веса t

	enc.printVector(u, "Информационный вектор u");
	enc.printVector(c, "Кодовый вектор");
	enc.printVector(y, format("Зашумленный кодовый вектор, ошибок - {}", t));

	vector<int> res = td.decode(y);						// результат декодирования
	vector<int> diff = enc.vectorDiff(c,res);			// проверка правильности декодирования

	enc.printVector(res, "Результат декодирования");
	enc.printVector(diff, "Разница между полученным вектором и исходным кодовым вектором");
}

// запуск iter_cnt итераций кодирования и декодирования случайного вектора и вычисление количества ошибочных декодирований
void testDecoderIter(int r1, int m1, int r2, int m2, int t, int iter_cnt) {

	ReedMallerTreeBuilder tb1(r1, m1);		// декодирующее дерево для кода Рида-Маллера 1
	ReedMallerTreeBuilder tb2(r2, m2);		// декодирующее дерево для кода Рида-Маллера 2
	TensorDecoder td(tb1.tree, tb2.tree);	// декодирующее дерево для тензорного произведения двух кодов Рида-Маллера

	cout << format("r1={}, m1={}, r2={}, m2={}, t={}", r1, m1, r2, m2, t) << endl;

	TensorEncoder enc;
	vector<vector<int>> G1 = enc.buildRMGeneratorMatrix(r1, m1);		// порождающая матрица для кода Рида-Маллера 1
	vector<vector<int>> G2 = enc.buildRMGeneratorMatrix(r2, m2);		// порождающая матрица для кода Рида-Маллера 2
	vector<vector<int>> Gt = enc.buildTensorMatrix(G1, G2);			// порождающая матрица для тензорного произведения двух кодов Рида-Маллера

	int k = (int)Gt.size(); // длина информационного вектора

	int s = 0; // количество ошибочных декодирований
	for (int i = 0; i < iter_cnt; ++i) {
       													
		vector<int> u = enc.generateRandomInfoVector(k);	// случайный информационный вектор
		vector<int> c = enc.encodeTensor(u, Gt);     		// кодовое слово
		vector<int> y = enc.addNoise(c, t);          		// добавили ошибку веса t

		vector<int> res = td.decode(y);						// результат декодирования
		vector<int> diff = enc.vectorDiff(c,res);			// проверка правильности декодирования

		bool error = false;
		for (int j=0; j<diff.size(); ++j)
			if (diff[j] != 0)
				error=true;

		if (error)
			s++;
	}

	cout << format("Ошибочных декодирований - {} из {}", s, iter_cnt) << endl << endl;
}

// тест декодера для разных кодов Рида-Маллера
void runDecoderTests(){
	
	testDecoderIter(0, 2, 0, 2, 7, 100); 
	testDecoderIter(0, 2, 0, 2, 1, 100); 
	testDecoderIter(0, 2, 1, 2, 3, 100); 
	testDecoderIter(0, 2, 1, 2, 1, 100); 
	
	testDecoderIter(1, 2, 0, 2, 3, 100); 
	testDecoderIter(1, 2, 0, 2, 1, 100); 
	testDecoderIter(1, 2, 1, 2, 1, 100); 

	testDecoderIter(0, 2, 0, 3, 15, 100); 
	testDecoderIter(0, 2, 1, 3, 7, 100);
	testDecoderIter(0, 2, 2, 3, 3, 100); 
	
	testDecoderIter(1, 2, 0, 3, 7, 100); 
	testDecoderIter(1, 2, 1, 3, 3, 100);
	testDecoderIter(1, 2, 2, 3, 1, 100); 

	testDecoderIter(0, 3, 0, 3, 31, 100); 
	testDecoderIter(0, 3, 1, 3, 15, 100);
	testDecoderIter(0, 3, 2, 3, 7, 100); 
	
	testDecoderIter(1, 3, 0, 3, 15, 100);
	testDecoderIter(1, 3, 1, 3, 7, 100);
	testDecoderIter(2, 3, 0, 3, 7, 100); 

}


int main()
{
	SetConsoleOutputCP(65001);
    SetConsoleCP(65001);

	int r1 = 1; int m1 = 2;
	int r2 = 1; int m2 = 3;
	int t = 3;
	
	// int iter_cnt = 100; 

	runDecoderRand(r1, m1, r2, m2, t);
	// testDecoderIter(r1, m1, r2, m2, t, iter_cnt);
	// runDecoderTests();

	system("pause");
	return 0;
}


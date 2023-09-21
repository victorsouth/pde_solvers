#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <clocale>
using namespace std;

class transport_equation_solver {
public:

	transport_equation_solver(int steps, double pipe_len, double step_len)
	{
		N = steps;
		n = static_cast<int>(pipe_len / step_len + 0.5) + 1;//считаем количество точек
	}


	vector<vector<double>> simple_moc_solver(int step, vector<double> layer_prev, vector<double> layer_curr, vector<double> oil_in, vector<double> oil_out, double flow)
	{
		int i;
		if (flow[i] > 0)
		{

			for (int l = 0; l < n - 1; l++)
				layer_curr[l + 1] = layer_prev[l];
			layer_curr[0] = oil_in[0][N - 1 - i];
		}
		else
		{
			for (int l = n - 1; l > 0; l--)
				layer_curr[l - 1] = layer_prev[l];
			layer_curr[n - i] = oil_out[0][i];

		}
	}

	int get_point_count(double pipe_len, double step_len) const {
		return n;
	}

	int get_steps() {
		return N;
	}

	vector<double> get_layer_prev() const {
		return layer_prev;
	}

	vector<double> get_layer_curr() const {
		return layer_curr;
	}

	void print_layers(int step, vector<double> layer_prev, vector<double> layer_curr) {
		ofstream fout("layers.txt");
		if (step == 0)
		{
			for (int j = 0; j < n; j++)
			{
				fout << layer_prev[j] << "\t";
			}
			fout << "\n";
		}
		for (int j = 0; j < n; j++)
		{
			fout << layer_curr[j] << "\t";
		}
		fout << "\n";
	}

private:
	int N;
	int n;

	vector<double> layer_prev;
	vector<double> layer_curr;
};
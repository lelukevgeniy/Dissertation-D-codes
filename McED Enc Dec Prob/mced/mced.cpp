// mced.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.


#include "../Codes/dcode.h"

int main(int argc, char* argv[])
{
	if (atoi(argv[1]) == 1)
	{
		//	argv[2] -- файл с параметрами кода
		//	argv[3] -- максимальное число итераций ISD 
		//	argv[4] -- число экспериментов
		DCode::test_EncryptionDecryption_ISD(argv[2], atoi(argv[3]), atoi(argv[4]));
	};

	if (atoi(argv[1]) == 2)
	{
		//	argv[2] -- файл с параметрами кода
		//	argv[3] -- число экспериментов
		DCode::test_EncryptionDecryption_MatrixDecoder(argv[2], atoi(argv[3]));
	};
}

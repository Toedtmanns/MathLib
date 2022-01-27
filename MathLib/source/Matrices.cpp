#include "../include/Maths.h"
#include <stdexcept>

namespace MathLib
{
	namespace Matrices
	{
		// MatrixI definition

		MatrixI::MatrixI(MatrixI&& other) noexcept
			: m_Rows(other.m_Rows), m_Columns(other.m_Columns), m_Matrix(other.m_Matrix)
		{
			other.m_Matrix = nullptr;
		}
		MatrixI::MatrixI(const MatrixI& other)
			: m_Rows(other.m_Rows), m_Columns(other.m_Columns), m_Matrix(new int* [other.m_Columns])
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				m_Matrix[col] = new int[m_Rows];
				memcpy(m_Matrix[col], other.m_Matrix[col], m_Rows * sizeof(int));
			}
		}
		MatrixI::MatrixI(const unsigned int& dim)
			: m_Rows(dim), m_Columns(dim), m_Matrix(new int*[dim])
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				m_Matrix[col] = new int[m_Rows]{0};
				m_Matrix[col][col] = 1;
			}
		}
		MatrixI::MatrixI(const unsigned int& columns, const unsigned int& rows)
			: m_Rows(rows), m_Columns(columns), m_Matrix(new int*[columns])
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				m_Matrix[col] = new int[m_Rows]{0};
			}
		}
		MatrixI::MatrixI(const unsigned int& columns, const unsigned int& rows, const int** const columnArray)
			: m_Rows(rows), m_Columns(columns), m_Matrix(new int*[columns])
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				m_Matrix[col] = new int[m_Rows];
				memcpy(m_Matrix[col], columnArray[col], m_Rows * sizeof(int));
			}
		}
		void MatrixI::SetNum(const unsigned int& column, const unsigned int& row, const int& value)
		{
			m_Matrix[column][row] = value;
		}
		void MatrixI::SetColumn(const unsigned int& column, const int* const content)
		{
			memcpy(m_Matrix[column], content, m_Rows * sizeof(int));
		}
		void MatrixI::SetMatrix(const int** const matrix)
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				memcpy(m_Matrix[col], matrix[col], m_Rows * sizeof(int));
			}
		}
		const int& MatrixI::GetNum(const unsigned int& column, const unsigned int& row) const
		{
			return m_Matrix[column][row];
		}
		int* MatrixI::GetColumn(const unsigned int& column) const
		{
			int* retCol = new int[m_Rows];
			memcpy(retCol, m_Matrix[column], m_Rows * sizeof(int));

			return retCol;
		}
		int** MatrixI::GetMatrix() const
		{
			int** retMat = new int* [m_Columns];
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				retMat[col] = new int[m_Rows];
				memcpy(retMat[col], m_Matrix[col], m_Rows * sizeof(int));
			}

			return retMat;
		}
		const unsigned int& MatrixI::GetRowCount() const
		{
			return m_Rows;
		}
		const unsigned int& MatrixI::GetColumnCount() const
		{
			return m_Columns;
		}
		void MatrixI::operator=(MatrixI&& other) noexcept
		{
			for (unsigned int col = 0; col < m_Columns; col++)
				delete[] m_Matrix[col];
			delete[] m_Matrix;

			m_Columns = other.m_Columns;
			m_Rows = other.m_Rows;
			m_Matrix = other.m_Matrix;
			other.m_Matrix = nullptr;
		}
		void MatrixI::operator=(const MatrixI& other)
		{
			m_Columns = other.m_Columns;
			m_Rows = other.m_Rows;

			for (unsigned int col = 0; col < other.m_Columns; col++)
			{
				delete[] m_Matrix[col];
				m_Matrix[col] = new int[m_Rows];
				memcpy(m_Matrix[col], other.m_Matrix[col], other.m_Rows * sizeof(int));
			}
		}
		MatrixI MatrixI::operator+(const int& value) const
		{
			MatrixI tempMat = MatrixI(m_Columns, m_Rows);
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] + value);
				}
			}

			return tempMat;
		}
		MatrixI MatrixI::operator-(const int& value) const
		{
			MatrixI tempMat = MatrixI(m_Columns, m_Rows);
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] + value);
				}
			}

			return tempMat;
		}
		MatrixI MatrixI::operator*(const int& value) const
		{
			MatrixI tempMat = MatrixI(m_Columns, m_Rows);
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] + value);
				}
			}

			return tempMat;
		}
		MatrixI MatrixI::operator/(const int& value) const
		{
			MatrixI tempMat = MatrixI(m_Columns, m_Rows);
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] + value);
				}
			}

			return tempMat;
		}
		MatrixI MatrixI::operator+(const MatrixI& other) const
		{
			if (m_Columns != other.GetColumnCount() || m_Rows != other.GetRowCount())
			{
				return MatrixI(0, 0);
			}
			MatrixI tempMat = MatrixI(m_Columns, m_Rows);

			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] + other.GetNum(col, row));
				}
			}

			return tempMat;
		}
		MatrixI MatrixI::operator-(const MatrixI& other) const
		{
			if (m_Columns != other.GetColumnCount() || m_Rows != other.GetRowCount())
			{
				return MatrixI(0, 0);
			}
			MatrixI tempMat = MatrixI(m_Columns, m_Rows);

			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] - other.GetNum(col, row));
				}
			}

			return tempMat;
		}
		MatrixI MatrixI::operator*(const MatrixI& other) const
		{
			if (m_Columns != other.GetRowCount())
				return MatrixI(0, 0);

			MatrixI tempMat = MatrixI(m_Rows, other.GetColumnCount());

			for (unsigned int col = 0; col < other.GetColumnCount(); col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					int num = 0;

					for (unsigned int calcCol = 0; calcCol < m_Columns; calcCol++)
					{
						num = num + m_Matrix[calcCol][row] * other.GetNum(calcCol, col);
					}

					tempMat.SetNum(row, col, num);
				}
			}

			return tempMat;
		}
		int* MatrixI::operator[](const unsigned int& column)
		{
			if (column < m_Columns)
				return m_Matrix[column];
		}
		MatrixI::~MatrixI()
		{
			if (m_Matrix == nullptr)
				return;
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				delete[] m_Matrix[col];
			}
			delete[] m_Matrix;
		}

		// MatrixF definition

		MatrixF::MatrixF(MatrixF&& other) noexcept
			: m_Rows(other.m_Rows), m_Columns(other.m_Columns), m_Matrix(other.m_Matrix)
		{
			other.m_Matrix = nullptr;
		}
		MatrixF::MatrixF(const MatrixF& other)
			: m_Rows(other.m_Rows), m_Columns(other.m_Columns), m_Matrix(new double* [other.m_Columns])
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				m_Matrix[col] = new double[m_Rows];
				memcpy(m_Matrix[col], other.m_Matrix[col], m_Rows * sizeof(double));
			}
		}
		MatrixF::MatrixF(const unsigned int& dim)
			: m_Rows(dim), m_Columns(dim), m_Matrix(new double* [dim])
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				m_Matrix[col] = new double[m_Rows] { 0 };
				m_Matrix[col][col] = 1;
			}
		}
		MatrixF::MatrixF(const unsigned int& columns, const unsigned int& rows)
			: m_Rows(rows), m_Columns(columns), m_Matrix(new double* [columns])
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				m_Matrix[col] = new double[m_Rows] { 0 };
			}
		}
		MatrixF::MatrixF(const unsigned int& columns, const unsigned int& rows, const double** const columnArray)
			: m_Rows(rows), m_Columns(columns), m_Matrix(new double* [columns])
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				m_Matrix[col] = new double[m_Rows];
				memcpy(m_Matrix[col], columnArray[col], m_Rows * sizeof(double));
			}
		}
		void MatrixF::SetNum(const unsigned int& column, const unsigned int& row, const double& value)
		{
			m_Matrix[column][row] = value;
		}
		void MatrixF::SetColumn(const unsigned int& column, const double* const content)
		{
			memcpy(m_Matrix[column], content, m_Rows * sizeof(double));
		}
		void MatrixF::SetMatrix(const double** const matrix)
		{
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				memcpy(m_Matrix[col], matrix[col], m_Rows * sizeof(double));
			}
		}
		const double& MatrixF::GetNum(const unsigned int& column, const unsigned int& row) const
		{
			return m_Matrix[column][row];
		}
		double* MatrixF::GetColumn(const unsigned int& column) const
		{
			double* retCol = new double[m_Rows];
			memcpy(retCol, m_Matrix[column], m_Rows * sizeof(double));

			return retCol;
		}
		double** MatrixF::GetMatrix() const
		{
			double** retMat = new double* [m_Columns];
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				retMat[col] = new double[m_Rows];
				memcpy(retMat[col], m_Matrix[col], m_Rows * sizeof(double));
			}

			return retMat;
		}
		const unsigned int& MatrixF::GetRowCount() const
		{
			return m_Rows;
		}
		const unsigned int& MatrixF::GetColumnCount() const
		{
			return m_Columns;
		}
		void MatrixF::operator=(MatrixF&& other) noexcept
		{
			for (unsigned int col = 0; col < m_Columns; col++)
				delete[] m_Matrix[col];
			delete[] m_Matrix;

			m_Columns = other.m_Columns;
			m_Rows = other.m_Rows;
			m_Matrix = other.m_Matrix;
			other.m_Matrix = nullptr;
		}
		void MatrixF::operator=(const MatrixF& other)
		{
			m_Columns = other.m_Columns;
			m_Rows = other.m_Rows;

			for (unsigned int col = 0; col < other.m_Columns; col++)
			{
				delete[] m_Matrix[col];
				memcpy(m_Matrix[col], other.m_Matrix[col], other.m_Rows * sizeof(double));
			}
		}
		MatrixF MatrixF::operator+(const double& value) const
		{
			MatrixF tempMat = MatrixF(m_Columns, m_Rows);
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] + value);
				}
			}

			return tempMat;
		}
		MatrixF MatrixF::operator-(const double& value) const
		{
			MatrixF tempMat = MatrixF(m_Columns, m_Rows);
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] - value);
				}
			}

			return tempMat;
		}
		MatrixF MatrixF::operator*(const double& value) const
		{
			MatrixF tempMat = MatrixF(m_Columns, m_Rows);
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] * value);
				}
			}

			return tempMat;
		}
		MatrixF MatrixF::operator/(const double& value) const
		{
			MatrixF tempMat = MatrixF(m_Columns, m_Rows);
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] / value);
				}
			}

			return tempMat;
		}
		MatrixF MatrixF::operator+(const MatrixF& other) const
		{
			if (m_Columns != other.GetColumnCount() || m_Rows != other.GetRowCount())
			{
				return MatrixF(0, 0);
			}
			MatrixF tempMat = MatrixF(m_Columns, m_Rows);

			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] + other.GetNum(col, row));
				}
			}

			return tempMat;
		}
		MatrixF MatrixF::operator-(const MatrixF& other) const
		{
			if (m_Columns != other.GetColumnCount() || m_Rows != other.GetRowCount())
			{
				return MatrixF(0, 0);
			}
			MatrixF tempMat = MatrixF(m_Columns, m_Rows);

			for (unsigned int col = 0; col < m_Columns; col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					tempMat.SetNum(col, row, m_Matrix[col][row] - other.GetNum(col, row));
				}
			}

			return tempMat;
		}
		MatrixF MatrixF::operator*(const MatrixF& other) const
		{
			if (m_Columns != other.GetRowCount())
				return MatrixF(0, 0);

			MatrixF tempMat = MatrixF(m_Rows, other.GetColumnCount());

			for (unsigned int col = 0; col < other.GetColumnCount(); col++)
			{
				for (unsigned int row = 0; row < m_Rows; row++)
				{
					double num = 0;

					for (unsigned int calcCol = 0; calcCol < m_Columns; calcCol++)
					{
						num = num + m_Matrix[calcCol][row] * other.GetNum(calcCol, col);
					}

					tempMat.SetNum(row, col, num);
				}
			}

			return tempMat;
		}
		double* MatrixF::operator[](const unsigned int& column)
		{
			if (column < m_Columns)
				return m_Matrix[column];
		}
		MatrixF::~MatrixF()
		{
			if (m_Matrix == nullptr)
				return;
			for (unsigned int col = 0; col < m_Columns; col++)
			{
				delete[] m_Matrix[col];
			}
			delete[] m_Matrix;
		}

		// General MatrixI functions

		void PrintContent(const MatrixI& mat)
		{
			for (unsigned int row = 0; row < mat.GetRowCount(); row++)
			{
				for (unsigned int col = 0; col < mat.GetColumnCount(); col++)
				{
					printf("%i\t", mat.GetNum(col, row));
				}
				printf("\n");
			}
		}
		void PrintProperties(const MatrixI& mat)
		{
			printf("Rows: %i, Columns: %i\n", mat.GetRowCount(), mat.GetColumnCount());
		}
		int MatrixGetDet(const MatrixI& mat)
		{
			if (mat.GetRowCount() != mat.GetColumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			int num = 0;
			unsigned int dim = mat.GetRowCount();

			if (dim == 2)
			{
				num = mat.GetNum(0, 0) * mat.GetNum(1, 1) - mat.GetNum(0, 1) * mat.GetNum(1, 0);
			}
			else if (dim == 3)
			{
				int pos1 = mat.GetNum(0, 0) * mat.GetNum(1, 1) * mat.GetNum(2, 2);
				int pos2 = mat.GetNum(1, 0) * mat.GetNum(2, 1) * mat.GetNum(0, 2);
				int pos3 = mat.GetNum(2, 0) * mat.GetNum(0, 1) * mat.GetNum(1, 2);

				int neg1 = mat.GetNum(0, 2) * mat.GetNum(1, 1) * mat.GetNum(2, 0);
				int neg2 = mat.GetNum(1, 2) * mat.GetNum(2, 1) * mat.GetNum(0, 0);
				int neg3 = mat.GetNum(2, 2) * mat.GetNum(0, 1) * mat.GetNum(1, 0);

				num = pos1 + pos2 + pos3 - neg1 - neg2 - neg3;
			}
			else
			{
				for (unsigned int mCol = 0; mCol < dim; mCol++)
				{
					MatrixI tempMat = MatrixI(dim - 1, dim - 1);
					unsigned int wCol = 0;

					for (unsigned int col = 0; col < dim; col++)
					{
						if (col == mCol)
						{
							continue;
						}
						for (unsigned int row = 0; row < dim - 1; row++)
						{
							tempMat.SetNum(row, wCol, mat.GetNum(col, row + 1));
						}
						wCol++;
					}

					if (mCol % 2 == 0)
						num += mat.GetNum(mCol, 0) * MatrixGetDet(tempMat);
					else
						num -= mat.GetNum(mCol, 0) * MatrixGetDet(tempMat);
				}
			}

			return num;
		}
		MatrixI MatrixOfMinors(const MatrixI& mat)
		{
			if (mat.GetRowCount() != mat.GetColumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			unsigned int dim = mat.GetRowCount();

			MatrixI resMat = MatrixI(dim, dim);

			for (unsigned int mRow = 0; mRow < dim; mRow++)
			{
				for (unsigned int mCol = 0; mCol < dim; mCol++)
				{
					MatrixI tempMat = MatrixI(dim - 1, dim - 1);

					unsigned int wRow = 0;

					for (unsigned int row = 0; row < dim; row++)
					{
						unsigned int wCol = 0;

						if (row == mRow)
							continue;

						for (unsigned int col = 0; col < dim; col++)
						{
							if (col == mCol)
								continue;

							tempMat.SetNum(wCol, wRow, mat.GetNum(col, row));
							wCol++;
						}
						wRow++;
					}

					resMat.SetNum(mCol, mRow, MatrixGetDet(tempMat));
				}
			}

			return resMat;
		}
		MatrixI MatrixOfCofactors(const MatrixI& mat)
		{
			unsigned int counter = 0;
			MatrixI resMat = MatrixI(mat.GetColumnCount(), mat.GetRowCount());

			for (unsigned int row = 0; row < mat.GetRowCount(); row++)
			{
				for (unsigned int col = 0; col < mat.GetColumnCount(); col++)
				{
					if (counter % 2 == 0)
						resMat.SetNum(col, row, mat.GetNum(col, row));
					else
						resMat.SetNum(col, row, -mat.GetNum(col, row));

					counter++;
				}
				if (mat.GetColumnCount() % 2 == 0)
					counter++;
			}

			return resMat;
		}
		MatrixI MatrixTranspose(const MatrixI& mat)
		{
			if (mat.GetRowCount() != mat.GetColumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			unsigned int dim = mat.GetRowCount();
			MatrixI resMat = MatrixI(dim, dim);

			for (unsigned int row = 0; row < dim; row++)
			{
				for (unsigned int col = 0; col < dim; col++)
				{
					resMat.SetNum(row, col, mat.GetNum(row, col));
				}
			}

			return resMat;
		}
		MatrixI MatrixAdjugate(const MatrixI& mat)
		{
			MatrixI retMat = MatrixOfCofactors(mat);
			return MatrixTranspose(retMat);
		}
		MatrixF MatrixInverse(const MatrixI& mat)
		{
			int determinant = MatrixGetDet(mat);

			unsigned int dim = mat.GetRowCount();
			MatrixF resMat = MatrixF(dim, dim);
			MatrixI tempMat = MatrixI(dim, dim);

			if (dim == 2)
			{
				resMat.SetNum(0, 0, mat.GetNum(1, 1));
				resMat.SetNum(1, 0, -1 * (double)mat.GetNum(1, 0));
				resMat.SetNum(0, 1, -1 * (double)mat.GetNum(0, 1));
				resMat.SetNum(1, 1, mat.GetNum(0, 0));

				resMat = resMat * (1 / MatrixGetDet(mat));
			}
			else
			{
				tempMat = MatrixOfMinors(mat);
				tempMat = MatrixOfCofactors(tempMat);
				tempMat = MatrixTranspose(tempMat);
				resMat = MatrixI2F(tempMat);
				resMat = resMat * (1 / determinant);
			}

			return resMat;
		}

		// General MatrixF functions

		void PrintContent(const MatrixF& mat)
		{
			for (unsigned int row = 0; row < mat.GetRowCount(); row++)
			{
				for (unsigned int col = 0; col < mat.GetColumnCount(); col++)
				{
					printf("%.3f\t", mat.GetNum(col, row));
				}
				printf("\n");
			}
		}
		void PrintProperties(const MatrixF& mat)
		{
			printf("Rows: %i, Columns: %i\n", mat.GetRowCount(), mat.GetColumnCount());
		}
		double MatrixGetDet(const MatrixF& mat)
		{
			if (mat.GetRowCount() != mat.GetColumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			double num = 0;
			unsigned int dim = mat.GetRowCount();

			if (dim == 2)
			{
				num = mat.GetNum(0, 0) * mat.GetNum(1, 1) - mat.GetNum(0, 1) * mat.GetNum(1, 0);
			}
			else if (dim == 3)
			{
				double pos1 = mat.GetNum(0, 0) * mat.GetNum(1, 1) * mat.GetNum(2, 2);
				double pos2 = mat.GetNum(1, 0) * mat.GetNum(2, 1) * mat.GetNum(0, 2);
				double pos3 = mat.GetNum(2, 0) * mat.GetNum(0, 1) * mat.GetNum(1, 2);

				double neg1 = mat.GetNum(0, 2) * mat.GetNum(1, 1) * mat.GetNum(2, 0);
				double neg2 = mat.GetNum(1, 2) * mat.GetNum(2, 1) * mat.GetNum(0, 0);
				double neg3 = mat.GetNum(2, 2) * mat.GetNum(0, 1) * mat.GetNum(1, 0);

				num = pos1 + pos2 + pos3 - neg1 - neg2 - neg3;
			}
			else
			{
				for (unsigned int mCol = 0; mCol < dim; mCol++)
				{
					MatrixF tempMat = MatrixF(dim - 1, dim - 1);
					unsigned int wCol = 0;

					for (unsigned int col = 0; col < dim; col++)
					{
						if (col == mCol)
						{
							continue;
						}
						for (unsigned int row = 0; row < dim - 1; row++)
						{
							tempMat.SetNum(row, wCol, mat.GetNum(col, row + 1));
						}
						wCol++;
					}

					if (mCol % 2 == 0)
						num += mat.GetNum(mCol, 0) * MatrixGetDet(tempMat);
					else
						num -= mat.GetNum(mCol, 0) * MatrixGetDet(tempMat);
				}
			}

			return num;
		}
		MatrixF MatrixOfMinors(const MatrixF& mat)
		{
			if (mat.GetRowCount() != mat.GetColumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			unsigned int dim = mat.GetRowCount();

			MatrixF resMat = MatrixF(dim, dim);

			for (unsigned int mRow = 0; mRow < dim; mRow++)
			{
				for (unsigned int mCol = 0; mCol < dim; mCol++)
				{
					MatrixF tempMat = MatrixF(dim - 1, dim - 1);

					unsigned int wRow = 0;

					for (unsigned int row = 0; row < dim; row++)
					{
						unsigned int wCol = 0;

						if (row == mRow)
							continue;

						for (unsigned int col = 0; col < dim; col++)
						{
							if (col == mCol)
								continue;

							tempMat.SetNum(wCol, wRow, mat.GetNum(col, row));
							wCol++;
						}
						wRow++;
					}

					resMat.SetNum(mCol, mRow, MatrixGetDet(tempMat));
				}
			}

			return resMat;
		}
		MatrixF MatrixOfCofactors(const MatrixF& mat)
		{
			unsigned int counter = 0;
			MatrixF resMat = MatrixF(mat.GetColumnCount(), mat.GetRowCount());

			for (unsigned int row = 0; row < mat.GetRowCount(); row++)
			{
				for (unsigned int col = 0; col < mat.GetColumnCount(); col++)
				{
					if (counter % 2 == 0)
						resMat.SetNum(col, row, mat.GetNum(col, row));
					else
						resMat.SetNum(col, row, -mat.GetNum(col, row));

					counter++;
				}
				if (mat.GetColumnCount() % 2 == 0)
					counter++;
			}

			return resMat;
		}
		MatrixF MatrixTranspose(const MatrixF& mat)
		{
			if (mat.GetRowCount() != mat.GetColumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			unsigned int dim = mat.GetRowCount();
			MatrixF resMat = MatrixF(dim, dim);

			for (unsigned int row = 0; row < dim; row++)
			{
				for (unsigned int col = 0; col < dim; col++)
				{
					resMat.SetNum(row, col, mat.GetNum(row, col));
				}
			}

			return resMat;
		}
		MatrixF MatrixAdjugate(const MatrixF& mat)
		{
			MatrixF retMat = MatrixOfCofactors(mat);
			return MatrixTranspose(retMat);
		}
		MatrixF MatrixInverse(const MatrixF& mat)
		{
			double determinant = MatrixGetDet(mat);

			unsigned int dim = mat.GetRowCount();
			MatrixF resMat = MatrixF(dim, dim);

			if (dim == 2)
			{
				resMat.SetNum(0, 0, mat.GetNum(1, 1));
				resMat.SetNum(1, 0, -1 * mat.GetNum(1, 0));
				resMat.SetNum(0, 1, -1 * mat.GetNum(0, 1));
				resMat.SetNum(1, 1, mat.GetNum(0, 0));

				resMat = resMat * (1 / MatrixGetDet(mat));
			}
			else
			{
				resMat = MatrixOfMinors(mat);
				resMat = MatrixOfCofactors(resMat);
				resMat = MatrixTranspose(resMat);
				resMat = resMat * (1 / determinant);
			}

			return resMat;
		}

		// Quaternion operations

		MatrixF MatrixRotate(const MatrixF& mat, const Complex::Quaternion& quat)
		{
			if (!MatrixIsSquare(mat, 4))
				throw std::invalid_argument("Matrix is not 4x4!");

			Complex::Quaternion xQuat(Primitives::Float3(mat.GetNum(0, 0), mat.GetNum(0, 1), mat.GetNum(0, 2)));
			Complex::Quaternion yQuat(Primitives::Float3(mat.GetNum(1, 0), mat.GetNum(1, 1), mat.GetNum(1, 2)));
			Complex::Quaternion zQuat(Primitives::Float3(mat.GetNum(2, 0), mat.GetNum(2, 1), mat.GetNum(2, 2)));

			xQuat = const_cast<Complex::Quaternion&>(quat).RotateQuaternion(xQuat);
			yQuat = const_cast<Complex::Quaternion&>(quat).RotateQuaternion(yQuat);
			zQuat = const_cast<Complex::Quaternion&>(quat).RotateQuaternion(zQuat);

			double col0[4] = {xQuat.i.num, xQuat.j.num, xQuat.k.num, mat.GetNum(0, 4)};
			double col1[4] = {yQuat.i.num, yQuat.j.num, yQuat.k.num, mat.GetNum(1, 4)};
			double col2[4] = {zQuat.i.num, zQuat.j.num, zQuat.k.num, mat.GetNum(2, 4)};

			MatrixF retMat = MatrixF(4);

			retMat.SetColumn(0, col0);
			retMat.SetColumn(1, col1);
			retMat.SetColumn(2, col2);
			retMat.SetColumn(3, mat.GetColumn(3));

			return retMat;
		}

		// Conversion and checking functions

		MatrixI MatrixF2I(const MatrixF& mat)
		{
			MatrixI resMat = MatrixI(mat.GetRowCount(), mat.GetColumnCount());

			for (unsigned int row = 0; row < resMat.GetRowCount(); row++)
			{
				for (unsigned int col = 0; col < resMat.GetColumnCount(); col++)
				{
					resMat.SetNum(row, col, mat.GetNum(col, row));
				}
			}

			return resMat;
		}
		MatrixF MatrixI2F(const MatrixI& mat)
		{
			MatrixF resMat = MatrixF(mat.GetRowCount(), mat.GetColumnCount());

			for (unsigned int row = 0; row < resMat.GetRowCount(); row++)
			{
				for (unsigned int col = 0; col < resMat.GetColumnCount(); col++)
				{
					resMat.SetNum(row, col, mat.GetNum(col, row));
				}
			}

			return resMat;
		}
		bool MatrixIsSquare(const MatrixI& mat, const int& dimension)
		{
			if (mat.GetRowCount() == mat.GetColumnCount())
			{
				if (dimension == -1)
				{
					return true;
				}
				else
				{
					return mat.GetRowCount() == dimension;
				}
			}
			return false;
		}
		bool MatrixIsSquare(const MatrixF& mat, const int& dimension)
		{
			if (mat.GetRowCount() == mat.GetColumnCount())
			{
				if (dimension == -1)
				{
					return true;
				}
				else
				{
					return mat.GetRowCount() == dimension;
				}
			}
			return false;
		}

		// Helper functions

		MatrixI Point2Matrix(const Primitives::Int2& point)
		{
			MatrixI retMat = MatrixI(2, 1);
			retMat.SetNum(0, 0, point.x);
			retMat.SetNum(1, 0, point.y);
			return retMat;
		}
		MatrixF Point2Matrix(const Primitives::Float2& point)
		{
			MatrixF retMat = MatrixF(2, 1);
			retMat.SetNum(0, 0, point.x);
			retMat.SetNum(1, 0, point.y);
			return retMat;
		}
		MatrixI Point2Matrix(const Primitives::Int3& point)
		{
			MatrixI retMat = MatrixI(3, 1);
			retMat.SetNum(0, 0, point.x);
			retMat.SetNum(1, 0, point.y);
			retMat.SetNum(2, 0, point.z);
			return retMat;
		}
		MatrixF Point2Matrix(const Primitives::Float3& point)
		{
			MatrixF retMat = MatrixF(3, 1);
			retMat.SetNum(0, 0, point.x);
			retMat.SetNum(1, 0, point.y);
			retMat.SetNum(2, 0, point.z);
			return retMat;
		}
		MatrixI TransformI2x2(const Primitives::Int2& p1, const Primitives::Int2& p2)
		{
			MatrixI retMat = MatrixI(2);
			retMat.SetNum(0, 0, p1.x);
			retMat.SetNum(0, 1, p1.y);
			retMat.SetNum(1, 0, p2.x);
			retMat.SetNum(1, 1, p2.y);

			return retMat;
		}
		MatrixF TransformF2x2(const Primitives::Float2& p1, const Primitives::Float2& p2)
		{
			MatrixF retMat = MatrixF(2);
			retMat.SetNum(0, 0, p1.x);
			retMat.SetNum(0, 1, p1.y);
			retMat.SetNum(1, 0, p2.x);
			retMat.SetNum(1, 1, p2.y);
			
			return retMat;
		}
	}
}
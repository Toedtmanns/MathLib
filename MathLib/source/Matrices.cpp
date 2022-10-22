#include "../include/MathLib/MathLib.hpp"
#include <stdexcept>

namespace MathLib
{
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
	double* MatrixF::GetArray() const
	{
		double* retArray = new double[m_Rows * m_Columns];

		double* rArrayInsert = retArray;
		for (unsigned int c = 0; c < m_Columns; c++)
		{
			memcpy(rArrayInsert, m_Matrix[c], m_Rows * sizeof(double));
			rArrayInsert += m_Rows;
		}
		return retArray;
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
		return m_Matrix[column];
	}
	const double* MatrixF::operator[](const unsigned int& column) const
	{
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
				resMat.SetNum(row, col, mat.GetNum(col, row));
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

	// Conversion and checking functions

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

	MatrixF Point2Matrix(const Float2& point)
	{
		MatrixF retMat = MatrixF(2, 1);
		retMat.SetNum(0, 0, point.x);
		retMat.SetNum(1, 0, point.y);
		return retMat;
	}
	MatrixF Point2Matrix(const Float3& point)
	{
		MatrixF retMat = MatrixF(3, 1);
		retMat.SetNum(0, 0, point.x);
		retMat.SetNum(1, 0, point.y);
		retMat.SetNum(2, 0, point.z);
		return retMat;
	}
	MatrixF TransformF2x2(const Float2& p1, const Float2& p2)
	{
		MatrixF retMat = MatrixF(2);
		retMat.SetNum(0, 0, p1.x);
		retMat.SetNum(0, 1, p1.y);
		retMat.SetNum(1, 0, p2.x);
		retMat.SetNum(1, 1, p2.y);
		
		return retMat;
	}
}
#include "../include/Maths.h"

namespace MathLib
{
	namespace Matrices
	{
		// RowI definition

		RowI::RowI(int length = 1, int val = 0)
			: length(length)
		{
			row = std::vector<int>(length);
			std::fill(row.begin(), row.end(), val);
		}
		RowI::RowI(std::vector<int> row)
			: length(row.size())
		{
			this->row = row;
		}
		int RowI::GetLength() const
		{
			return length;
		}
		int RowI::GetAt(int index)
		{
			return row[index];
		}
		std::vector<int> RowI::GetRow() const
		{
			std::vector<int> ret = std::vector<int>();
			for (int i = 0; i < length; i++)
			{
				ret.push_back(row[i]);
			}

			return ret;
		}
		void RowI::SetRow(std::vector<int> *content)
		{
			int cpyLen = content->size();
			if (cpyLen > length)
				cpyLen = length;

			for (int i = 0; i < cpyLen - 1; i++)
				row[i] = content->at(i);
		}
		void RowI::SetNum(int collumn, int value)
		{
			row[collumn] = value;
		}

		// RowF definition

		RowF::RowF(int length = 1, double val = 0)
			: length(length)
		{
			row = std::vector<double>(length);
			std::fill(row.begin(), row.end(), val);
		}
		RowF::RowF(std::vector<double> row)
			: length(row.size())
		{
			this->row = row;
		}
		int RowF::GetLength() const
		{
			return length;
		}
		double RowF::GetAt(int index)
		{
			return row[index];
		}
		std::vector<int> RowF::GetRow() const
		{
			std::vector<int> ret = std::vector<int>();
			for (int i = 0; i < length; i++)
			{
				ret.push_back(row[i]);
			}

			return ret;
		}
		void RowF::SetRow(std::vector<double> *content)
		{
			int cpyLen = content->size();
			if (cpyLen > length)
				cpyLen = length;

			for (int i = 0; i < cpyLen - 1; i++)
				row[i] = content->at(i);
		}
		void RowF::SetNum(int collumn, double value)
		{
			row[collumn] = value;
		}

		// MatrixI definition

		MatrixI::MatrixI(int rows = 2, int collumns = 2)
			: rows(rows), collumns(collumns)
		{
			matrix = std::vector<RowI>(rows);
			for (int i = 0; i < rows; i++)
			{
				RowI temp = RowI(collumns);
				matrix[i] = temp;
			}
		}
		MatrixI::MatrixI(std::vector<RowI> content)
			: rows(content.size()), collumns(content[0].GetLength())
		{
			matrix = content;
		}
		int MatrixI::GetRowCount() const
		{
			return rows;
		}
		int MatrixI::GetCollumnCount() const
		{
			return collumns;
		}
		RowI MatrixI::GetRow(int index)
		{
			return matrix[index];
		}
		void MatrixI::SetRow(int index, RowI row)
		{
			matrix[index] = row;
		}
		void MatrixI::SetNum(int row, int collumn, int value)
		{
			matrix[row].SetNum(collumn, value);
		}

		// MatrixF definition

		MatrixF::MatrixF(int dim)
			: rows(dim), collumns(dim)
		{
			matrix = std::vector<RowF>(dim);
			for (int i = 0; i < dim; i++)
			{
				RowF temp = RowF(dim);
				temp.SetNum(i, 1);
				matrix[i] = temp;
			}
		}
		MatrixF::MatrixF(int rows = 2, int collumns = 2)
			: rows(rows), collumns(collumns)
		{
			matrix = std::vector<RowF>(rows);
			for (int i = 0; i < rows; i++)
			{
				RowF temp = RowF(collumns);
				matrix[i] = temp;
			}
		}
		MatrixF::MatrixF(std::vector<RowF> content)
			: rows(content.size()), collumns(content[0].GetLength())
		{
			matrix = content;
		}
		int MatrixF::GetRowCount() const
		{
			return rows;
		}
		int MatrixF::GetCollumnCount() const
		{
			return collumns;
		}
		RowF MatrixF::GetRow(int index)
		{
			return matrix[index];
		}
		void MatrixF::SetRow(int index, RowF row)
		{
			matrix[index] = row;
		}
		void MatrixF::SetNum(int row, int collumn, double value)
		{
			matrix[row].SetNum(collumn, value);
		}
		MatrixF::~MatrixF()
		{
			matrix.clear();
		}

		// General MatrixI functions

		void PrintContent(MatrixI *mat)
		{
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					printf("%i\t", mat->GetRow(row).GetAt(col));
				}
				printf("\n");
			}
		}
		void PrintProperties(MatrixI *mat)
		{
			printf("Rows: %i, Collumns: %i\n", mat->GetRowCount(), mat->GetCollumnCount());
		}
		int MatrixGetDet(MatrixI *mat)
		{
			if (mat->GetRowCount() != mat->GetCollumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			int num = 0;
			int dim = mat->GetRowCount();

			if (dim == 2)
			{
				num = mat->GetRow(0).GetAt(0) * mat->GetRow(1).GetAt(1) - mat->GetRow(1).GetAt(0) * mat->GetRow(0).GetAt(1);
			}
			else if (dim == 3)
			{
				int pos1 = mat->GetRow(0).GetAt(0) * mat->GetRow(1).GetAt(1) * mat->GetRow(2).GetAt(2);
				int pos2 = mat->GetRow(0).GetAt(1) * mat->GetRow(1).GetAt(2) * mat->GetRow(2).GetAt(0);
				int pos3 = mat->GetRow(0).GetAt(2) * mat->GetRow(1).GetAt(0) * mat->GetRow(2).GetAt(1);

				int neg1 = mat->GetRow(2).GetAt(0) * mat->GetRow(1).GetAt(1) * mat->GetRow(0).GetAt(2);
				int neg2 = mat->GetRow(2).GetAt(1) * mat->GetRow(1).GetAt(2) * mat->GetRow(0).GetAt(0);
				int neg3 = mat->GetRow(2).GetAt(2) * mat->GetRow(1).GetAt(0) * mat->GetRow(0).GetAt(1);

				num = pos1 + pos2 + pos3 - neg1 - neg2 - neg3;
				//num = pos1 - neg2 + pos2 - neg3 + pos3 - neg1;
			}
			else
			{
				for (int mCol = 0; mCol < dim; mCol++)
				{
					MatrixI *tempMat = new MatrixI(dim - 1, dim - 1);
					int wCol = 0;

					for (int col = 0; col < dim; col++)
					{
						if (col == mCol)
						{
							continue;
						}
						for (int row = 0; row < dim - 1; row++)
						{
							tempMat->SetNum(row, wCol, mat->GetRow(row + 1).GetAt(col));
						}
						wCol++;
					}

					if (mCol % 2 == 0)
						num += mat->GetRow(0).GetAt(mCol) * MatrixGetDet(tempMat);
					else
						num -= mat->GetRow(0).GetAt(mCol) * MatrixGetDet(tempMat);
				}
			}

			return num;
		}
		MatrixI MatrixAdd(MatrixI *mat, int value)
		{
			MatrixI *tempMat = new MatrixI(mat->GetRowCount(), mat->GetCollumnCount());
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				RowI tempRow = mat->GetRow(row);
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					tempRow.SetNum(col, tempRow.GetAt(col) + value);
				}
				tempMat->SetRow(row, tempRow);
			}

			return *tempMat;
		}
		MatrixI MatrixSub(MatrixI *mat, int value)
		{
			MatrixI *tempMat = new MatrixI(mat->GetRowCount(), mat->GetCollumnCount());
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				RowI tempRow = mat->GetRow(row);
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					tempRow.SetNum(col, tempRow.GetAt(col) - value);
				}
				tempMat->SetRow(row, tempRow);
			}

			return *tempMat;
		}
		MatrixI MatrixMult(MatrixI *mat, int value)
		{
			MatrixI *tempMat = new MatrixI(mat->GetRowCount(), mat->GetCollumnCount());
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				RowI tempRow = mat->GetRow(row);
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					tempRow.SetNum(col, tempRow.GetAt(col) * value);
				}
				tempMat->SetRow(row, tempRow);
			}

			return *tempMat;
		}
		MatrixI MatrixDiv(MatrixI *mat, int value)
		{
			MatrixI *tempMat = new MatrixI(mat->GetRowCount(), mat->GetCollumnCount());
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				RowI tempRow = mat->GetRow(row);
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					tempRow.SetNum(col, tempRow.GetAt(col) / value);
				}
				tempMat->SetRow(row, tempRow);
			}

			return *tempMat;
		}
		MatrixI MatrixAdd(MatrixI *mat1, MatrixI *mat2)
		{
			if (mat1->GetCollumnCount() != mat2->GetCollumnCount() || mat1->GetRowCount() != mat2->GetRowCount())
			{
				return MatrixI(0, 0);
			}
			MatrixI tempMat = MatrixI(mat1->GetRowCount(), mat1->GetCollumnCount());

			for (int row = 0; row < mat1->GetRowCount(); row++)
			{
				for (int col = 0; col < mat1->GetCollumnCount(); col++)
				{
					int num = mat1->GetRow(row).GetAt(col) + mat2->GetRow(row).GetAt(col);
					tempMat.SetNum(row, col, num);
				}
			}

			return tempMat;
		}
		MatrixI MatrixSub(MatrixI *mat1, MatrixI *mat2)
		{
			if (mat1->GetCollumnCount() != mat2->GetCollumnCount() || mat1->GetRowCount() != mat2->GetRowCount())
			{
				return MatrixI(0, 0);
			}
			MatrixI tempMat = MatrixI(mat1->GetRowCount(), mat1->GetCollumnCount());

			for (int row = 0; row < mat1->GetRowCount(); row++)
			{
				for (int col = 0; col < mat1->GetCollumnCount(); col++)
				{
					int num = mat1->GetRow(row).GetAt(col) - mat2->GetRow(row).GetAt(col);
					tempMat.SetNum(row, col, num);
				}
			}

			return tempMat;
		}
		MatrixI MatrixMult(MatrixI *mat1, MatrixI *mat2)
		{
			if (mat1->GetCollumnCount() != mat2->GetRowCount())
				return MatrixI(0, 0);

			MatrixI tempMat = MatrixI(mat1->GetRowCount(), mat2->GetCollumnCount());

			for (int row = 0; row < mat1->GetRowCount(); row++)
			{
				for (int col = 0; col < mat1->GetRowCount(); col++)
				{
					int num = 0;

					for (int calcCol = 0; calcCol < mat1->GetCollumnCount(); calcCol++)
					{
						num = num + mat1->GetRow(row).GetAt(calcCol) * mat2->GetRow(calcCol).GetAt(col);
					}

					tempMat.SetNum(row, col, num);
				}
			}

			return tempMat;
		}
		MatrixI MatrixOfMinors(MatrixI *mat)
		{
			if (mat->GetRowCount() != mat->GetCollumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			int dim = mat->GetRowCount();

			MatrixI resMat = MatrixI(dim, dim);

			for (int mRow = 0; mRow < dim; mRow++)
			{
				for (int mCol = 0; mCol < dim; mCol++)
				{
					MatrixI *tempMat = new MatrixI(dim - 1, dim - 1);

					int wRow = 0;

					for (int row = 0; row < dim; row++)
					{
						int wCol = 0;

						if (row == mRow)
							continue;

						for (int col = 0; col < dim; col++)
						{
							if (col == mCol)
								continue;

							tempMat->SetNum(wRow, wCol, mat->GetRow(row).GetAt(col));
							wCol++;
						}
						wRow++;
					}

					resMat.SetNum(mRow, mCol, MatrixGetDet(tempMat));
				}
			}

			return resMat;
		}
		MatrixI MatrixOfCofactors(MatrixI *mat)
		{
			int counter = 0;
			MatrixI resMat = MatrixI(mat->GetRowCount(), mat->GetCollumnCount());

			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					if (counter % 2 == 0)
						resMat.SetNum(row, col, mat->GetRow(row).GetAt(col));
					else
						resMat.SetNum(row, col, -mat->GetRow(row).GetAt(col));

					counter++;
				}
				if (mat->GetCollumnCount() % 2 == 0)
					counter++;
			}

			return resMat;
		}
		MatrixI MatrixAdjugate(MatrixI *mat)
		{
			if (mat->GetRowCount() != mat->GetCollumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			int dim = mat->GetRowCount();
			MatrixI resMat = MatrixI(dim, dim);

			for (int row = 0; row < dim; row++)
			{
				for (int col = 0; col < dim; col++)
				{
					resMat.SetNum(row, col, mat->GetRow(col).GetAt(row));
				}
			}

			return resMat;
		}
		MatrixF MatrixInverse(MatrixI *mat)
		{
			int determinant = MatrixGetDet(mat);

			int dim = mat->GetRowCount();
			MatrixF resMat = MatrixF(dim, dim);
			MatrixI tempMat = MatrixI(dim, dim);

			if (dim == 2)
			{
				resMat.SetNum(0, 0, mat->GetRow(1).GetAt(1));
				resMat.SetNum(0, 1, -1 * (double)mat->GetRow(0).GetAt(1));
				resMat.SetNum(1, 0, -1 * (double)mat->GetRow(1).GetAt(0));
				resMat.SetNum(1, 1, mat->GetRow(0).GetAt(0));

				MatrixF *pRMat = &resMat;
				resMat = MatrixMult(pRMat, 1 / MatrixGetDet(mat));
				pRMat = nullptr;
			}
			else
			{
				MatrixI *pTMat = &tempMat;
				MatrixF *pRMat = &resMat;
				tempMat = MatrixOfMinors(mat);
				tempMat = MatrixOfCofactors(pTMat);
				tempMat = MatrixAdjugate(pTMat);
				resMat = MatrixI2F(pTMat);
				resMat = MatrixMult(pRMat, 1 / determinant);
			}

			return resMat;
		}

		// General MatrixF functions

		void PrintContent(MatrixF *mat)
		{
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					printf("%.3f\t", mat->GetRow(row).GetAt(col));
				}
				printf("\n");
			}
		}
		void PrintProperties(MatrixF *mat)
		{
			printf("Rows: %i, Collumns: %i\n", mat->GetRowCount(), mat->GetCollumnCount());
		}
		double MatrixGetDet(MatrixF *mat)
		{
			if (mat->GetRowCount() != mat->GetCollumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			double num = 0;
			int dim = mat->GetRowCount();

			if (dim == 2)
			{
				num = mat->GetRow(0).GetAt(0) * mat->GetRow(1).GetAt(1) - mat->GetRow(1).GetAt(0) * mat->GetRow(0).GetAt(1);
			}
			else if (dim == 3)
			{
				double pos1 = mat->GetRow(0).GetAt(0) * mat->GetRow(1).GetAt(1) * mat->GetRow(2).GetAt(2);
				double pos2 = mat->GetRow(0).GetAt(1) * mat->GetRow(1).GetAt(2) * mat->GetRow(2).GetAt(0);
				double pos3 = mat->GetRow(0).GetAt(2) * mat->GetRow(1).GetAt(0) * mat->GetRow(2).GetAt(1);

				double neg1 = mat->GetRow(2).GetAt(0) * mat->GetRow(1).GetAt(1) * mat->GetRow(0).GetAt(2);
				double neg2 = mat->GetRow(2).GetAt(1) * mat->GetRow(1).GetAt(2) * mat->GetRow(0).GetAt(0);
				double neg3 = mat->GetRow(2).GetAt(2) * mat->GetRow(1).GetAt(0) * mat->GetRow(0).GetAt(1);

				num = pos1 + pos2 + pos3 - neg1 - neg2 - neg3;
			}
			else
			{
				for (int mCol = 0; mCol < dim; mCol++)
				{
					MatrixF *tempMat = new MatrixF(dim - 1, dim - 1);
					int wCol = 0;

					for (int col = 0; col < dim; col++)
					{
						if (col == mCol)
						{
							continue;
						}
						for (int row = 0; row < dim - 1; row++)
						{
							tempMat->SetNum(row, wCol, mat->GetRow(row + 1).GetAt(col));
						}
						wCol++;
					}

					if (mCol % 2 == 0)
						num += mat->GetRow(0).GetAt(mCol) * MatrixGetDet(tempMat);
					else
						num -= mat->GetRow(0).GetAt(mCol) * MatrixGetDet(tempMat);

					delete tempMat;
				}
			}

			return num;
		}
		MatrixF MatrixAdd(MatrixF *mat, double value)
		{
			MatrixF tempMat = MatrixF(mat->GetRowCount(), mat->GetCollumnCount());
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				RowF tempRow = mat->GetRow(row);
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					tempRow.SetNum(col, tempRow.GetAt(col) + value);
				}
				tempMat.SetRow(row, tempRow);
			}

			return tempMat;
		}
		MatrixF MatrixSub(MatrixF *mat, double value)
		{
			MatrixF tempMat = MatrixF(mat->GetRowCount(), mat->GetCollumnCount());
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				RowF tempRow = mat->GetRow(row);
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					tempRow.SetNum(col, tempRow.GetAt(col) - value);
				}
				tempMat.SetRow(row, tempRow);
			}

			return tempMat;
		}
		MatrixF MatrixMult(MatrixF *mat, double value)
		{
			MatrixF tempMat = MatrixF(mat->GetRowCount(), mat->GetCollumnCount());
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				RowF tempRow = mat->GetRow(row);
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					tempRow.SetNum(col, tempRow.GetAt(col) * value);
				}
				tempMat.SetRow(row, tempRow);
			}

			return tempMat;
		}
		MatrixF MatrixDiv(MatrixF *mat, double value)
		{
			MatrixF *tempMat = new MatrixF(mat->GetRowCount(), mat->GetCollumnCount());
			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				RowF tempRow = mat->GetRow(row);
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					tempRow.SetNum(col, tempRow.GetAt(col) / value);
				}
				tempMat->SetRow(row, tempRow);
			}

			return *tempMat;
		}
		MatrixF MatrixAdd(MatrixF *mat1, MatrixF *mat2)
		{
			if (mat1->GetCollumnCount() != mat2->GetCollumnCount() || mat1->GetRowCount() != mat2->GetRowCount())
			{
				return MatrixF(0, 0);
			}
			MatrixF tempMat = MatrixF(mat1->GetRowCount(), mat1->GetCollumnCount());

			for (int row = 0; row < mat1->GetRowCount(); row++)
			{
				for (int col = 0; col < mat1->GetCollumnCount(); col++)
				{
					double num = mat1->GetRow(row).GetAt(col) + mat2->GetRow(row).GetAt(col);
					tempMat.SetNum(row, col, num);
				}
			}

			return tempMat;
		}
		MatrixF MatrixSub(MatrixF *mat1, MatrixF *mat2)
		{
			if (mat1->GetCollumnCount() != mat2->GetCollumnCount() || mat1->GetRowCount() != mat2->GetRowCount())
			{
				return MatrixF(0, 0);
			}
			MatrixF tempMat = MatrixF(mat1->GetRowCount(), mat1->GetCollumnCount());

			for (int row = 0; row < mat1->GetRowCount(); row++)
			{
				for (int col = 0; col < mat1->GetCollumnCount(); col++)
				{
					double num = mat1->GetRow(row).GetAt(col) - mat2->GetRow(row).GetAt(col);
					tempMat.SetNum(row, col, num);
				}
			}

			return tempMat;
		}
		MatrixF MatrixMult(MatrixF *mat1, MatrixF *mat2)
		{
			if (mat1->GetCollumnCount() != mat2->GetRowCount())
				return MatrixF(0, 0);

			MatrixF tempMat = MatrixF(mat1->GetRowCount(), mat2->GetCollumnCount());

			for (int row = 0; row < mat1->GetRowCount(); row++)
			{
				for (int col = 0; col < mat1->GetRowCount(); col++)
				{
					double num = 0;

					for (int calcCol = 0; calcCol < mat1->GetCollumnCount(); calcCol++)
					{
						num = num + mat1->GetRow(row).GetAt(calcCol) * mat2->GetRow(calcCol).GetAt(col);
					}

					tempMat.SetNum(row, col, num);
				}
			}

			return tempMat;
		}
		MatrixF MatrixOfMinors(MatrixF *mat)
		{
			if (mat->GetRowCount() != mat->GetCollumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			int dim = mat->GetRowCount();

			MatrixF resMat = MatrixF(dim, dim);

			for (int mRow = 0; mRow < dim; mRow++)
			{
				for (int mCol = 0; mCol < dim; mCol++)
				{
					MatrixF *tempMat = new MatrixF(dim - 1, dim - 1);

					int wRow = 0;

					for (int row = 0; row < dim; row++)
					{
						int wCol = 0;

						if (row == mRow)
							continue;

						for (int col = 0; col < dim; col++)
						{
							if (col == mCol)
								continue;

							tempMat->SetNum(wRow, wCol, mat->GetRow(row).GetAt(col));
							wCol++;
						}
						wRow++;
					}

					resMat.SetNum(mRow, mCol, MatrixGetDet(tempMat));
					delete tempMat;
				}
			}

			return resMat;
		}
		MatrixF MatrixOfCofactors(MatrixF *mat)
		{
			int counter = 0;
			MatrixF resMat = MatrixF(mat->GetRowCount(), mat->GetCollumnCount());

			for (int row = 0; row < mat->GetRowCount(); row++)
			{
				for (int col = 0; col < mat->GetCollumnCount(); col++)
				{
					if (counter % 2 == 0)
						resMat.SetNum(row, col, mat->GetRow(row).GetAt(col));
					else
						resMat.SetNum(row, col, -mat->GetRow(row).GetAt(col));

					counter++;
				}
				if (mat->GetCollumnCount() % 2 == 0)
					counter++;
			}

			return resMat;
		}
		MatrixF MatrixAdjugate(MatrixF *mat)
		{
			if (mat->GetRowCount() != mat->GetCollumnCount())
				throw std::invalid_argument("Matrix isn't square!");

			int dim = mat->GetRowCount();
			MatrixF resMat = MatrixF(dim, dim);

			for (int row = 0; row < dim; row++)
			{
				for (int col = 0; col < dim; col++)
				{
					resMat.SetNum(row, col, mat->GetRow(col).GetAt(row));
				}
			}

			return resMat;
		}
		MatrixF MatrixInverse(MatrixF *mat)
		{
			double determinant = MatrixGetDet(mat);

			int dim = mat->GetRowCount();
			MatrixF resMat = MatrixF(dim, dim);

			if (dim == 2)
			{
				resMat.SetNum(0, 0, mat->GetRow(1).GetAt(1));
				resMat.SetNum(0, 1, -1 * (double)mat->GetRow(0).GetAt(1));
				resMat.SetNum(1, 0, -1 * (double)mat->GetRow(1).GetAt(0));
				resMat.SetNum(1, 1, mat->GetRow(0).GetAt(0));

				resMat = MatrixMult(&resMat, 1 / MatrixGetDet(mat));
			}
			else
			{
				resMat = MatrixOfMinors(mat);
				resMat = MatrixOfCofactors(&resMat);
				resMat = MatrixAdjugate(&resMat);
				resMat = MatrixMult(&resMat, 1 / determinant);
			}

			return resMat;
		}
		
		// Conversion and checking functions

		MatrixI MatrixF2I(MatrixF *mat)
		{
			MatrixI resMat = MatrixI(mat->GetRowCount(), mat->GetCollumnCount());

			for (int row = 0; row < resMat.GetRowCount(); row++)
			{
				for (int col = 0; col < resMat.GetCollumnCount(); col++)
				{
					resMat.SetNum(row, col, mat->GetRow(row).GetAt(col));
				}
			}

			return resMat;
		}
		MatrixF MatrixI2F(MatrixI *mat)
		{
			MatrixF resMat = MatrixF(mat->GetRowCount(), mat->GetCollumnCount());

			for (int row = 0; row < resMat.GetRowCount(); row++)
			{
				for (int col = 0; col < resMat.GetCollumnCount(); col++)
				{
					resMat.SetNum(row, col, mat->GetRow(row).GetAt(col));
				}
			}

			return resMat;
		}
		bool MatrixIsSquare(const MatrixI &mat, int dimension)
		{
			if (mat.GetRowCount() == mat.GetCollumnCount())
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
		bool MatrixIsSquare(const MatrixF &mat, int dimension)
		{
			if (mat.GetRowCount() == mat.GetCollumnCount())
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

		MatrixI *MatrixI2x2(int r1c1, int r1c2, int r2c1, int r2c2)
		{
			return new MatrixI(std::vector<RowI>{RowI(std::vector<int>{r1c1, r1c2}), RowI(std::vector<int>{r2c1, r2c2})});
		}
		MatrixI *TransformI2x2(Primitives::Point2D c1, Primitives::Point2D c2)
		{
			return new MatrixI(std::vector<RowI>{RowI(std::vector<int>{(int)c1.x, (int)c2.x}), RowI(std::vector<int>{(int)c1.y, (int)c2.y})});
		}
		MatrixI *TransformI2x2()
		{
			return new MatrixI(std::vector<RowI>{RowI(std::vector<int>{1, 0}), RowI(std::vector<int>{0, 1})});
		}
		MatrixF *MatrixF2x2(double r1c1, double r1c2, double r2c1, double r2c2)
		{
			return new MatrixF(std::vector<RowF>{RowF(std::vector<double>{r1c1, r1c2}), RowF(std::vector<double>{r2c1, r2c2})});
		}
		MatrixF *TransformF2x2()
		{
			return new MatrixF(std::vector<RowF>{RowF(std::vector<double>{1, 0}), RowF(std::vector<double>{0, 1})});
		}
		MatrixF *TransformF2x2(Primitives::Point2D c1, Primitives::Point2D c2)
		{
			return new MatrixF(std::vector<RowF>{RowF(std::vector<double>{c1.x, c2.x}), RowF(std::vector<double>{c1.y, c2.y})});
		}
	}
}
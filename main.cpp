#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <fstream>

class StarAlignment
{
private:
    std::vector<std::string> sequences;
    std::string star_sequence;
    std::vector<std::vector<int>> pair_alignment_scores;
    std::vector<std::vector<std::pair<std::string, std::string>>> pair_optimal_alignments;
    int star_sequence_pos;

    std::pair<int, std::pair<std::string, std::string>> get_pair_alignment_score_and_sequence(const std::string & s, const std::string & t)
    {
        int rows = s.size() + 1;
        int cols = t.size() + 1;

        std::vector<int> diag_vector = {-2, -2};
        std::vector<int> prev_diag_vector = {0};
        std::vector<int> new_diag_vector;

        std::vector<std::pair<std::string, std::string>> alignment_dv = {std::make_pair("_", "" + t[0]), std::make_pair("" + s[0], "_")};
        std::vector<std::pair<std::string, std::string>> alignment_pdv = {std::make_pair("", "")};
        std::vector<std::pair<std::string, std::string>> alignment_ndv;

        for (int diagonal = 2; diagonal < rows + cols - 1; diagonal++)
        {
            int start_row = std::max(0, diagonal - cols + 1);
            int start_col = std::min(diagonal, cols - 1);

            int diag_len = std::min(rows - start_row - 1, start_col) + 1;
             
            int pos_dif = 0;
            if (start_col == cols - 1 && start_row >= 2) pos_dif = 1;
            else if (start_col == cols - 1 && start_row == 1) pos_dif = 0; 
            else pos_dif = -1;

            int left_dif = 0;
            if (start_row > 0) left_dif = 1;

            int top_dif = -1;
            if (start_col == cols - 1 && start_row > 0) top_dif = 0;

            new_diag_vector.resize(diag_len);
            alignment_ndv.resize(diag_len);

            #pragma omp parallel for
            for (int k = 0; k < diag_len; k++)
            {
                int i = start_row +  k, j = start_col - k;
                if (i == 0)
                {
                    new_diag_vector[k] = diag_vector[k + left_dif] - 2;
                    alignment_ndv[k] = std::make_pair(alignment_dv[k + left_dif].first + '_', alignment_dv[k + left_dif].second + t[j - 1]);
                    continue;
                }
                if (j == 0)
                {
                    new_diag_vector[k] = diag_vector[k + top_dif] - 2;
                    alignment_ndv[k] = std::make_pair(alignment_dv[k + top_dif].first + s[i - 1], alignment_dv[k + top_dif].second + '_');
                    continue;
                }

                int new_val;
                if (s[i - 1] == t[j - 1])
                {
                    new_val = std::max({diag_vector[k + top_dif] - 2, diag_vector[k + left_dif] - 2, prev_diag_vector[k + pos_dif] + 1});
                }
                else
                {
                    new_val = std::max({diag_vector[k + top_dif] - 2, diag_vector[k + left_dif] - 2, prev_diag_vector[k + pos_dif] - 1});
                }

                if ((new_val == prev_diag_vector[k + pos_dif] + 1 && s[i - 1] == t[j - 1]) || (new_val  == prev_diag_vector[k + pos_dif] - 1 && s[i - 1] != t[j - 1]))
                {
                    alignment_ndv[k] = std::make_pair(alignment_pdv[k + pos_dif].first + s[i - 1], alignment_pdv[k + pos_dif].second + t[j - 1]);
                }
                else if (new_val == diag_vector[k + top_dif] - 2)
                {
                    alignment_ndv[k] = std::make_pair(alignment_dv[k + top_dif].first + s[i - 1], alignment_dv[k + top_dif].second + '_');
                }
                else
                {
                    alignment_ndv[k] = std::make_pair(alignment_dv[k + left_dif].first + '_', alignment_dv[k + left_dif].second + t[j - 1]);
                }

                new_diag_vector[k] = new_val;
            }
            prev_diag_vector = diag_vector;
            diag_vector = new_diag_vector;
            new_diag_vector.clear();


            alignment_pdv = alignment_dv;
            alignment_dv = alignment_ndv;
            alignment_ndv.clear();
        }

        return std::make_pair(diag_vector[0], alignment_dv[0]);
    }

    std::string find_star_sequence()
    {
        this->pair_alignment_scores = std::vector<std::vector<int>>(sequences.size(), std::vector<int>(sequences.size(), 0));
        this->pair_optimal_alignments = std::vector<std::vector<std::pair<std::string, std::string>>> (sequences.size(), std::vector<std::pair<std::string, std::string>>(sequences.size()));

        std::vector<int> score_sums(sequences.size(), 0);
        for (int i = 0; i < sequences.size(); i++)
        {
            for (int j = i + 1; j < sequences.size(); j++)
            {
                auto score_and_alignment = get_pair_alignment_score_and_sequence(sequences[i], sequences[j]);
                this->pair_alignment_scores[i][j] = score_and_alignment.first;
                this->pair_alignment_scores[j][i] = this->pair_alignment_scores[i][j];
                score_sums[i] += this->pair_alignment_scores[i][j];
                score_sums[j] += this->pair_alignment_scores[i][j];

                this->pair_optimal_alignments[i][j] = score_and_alignment.second;
                this->pair_optimal_alignments[j][i] = std::make_pair(score_and_alignment.second.second, score_and_alignment.second.first);
            }
        }


        auto star_pos = std::max_element(score_sums.begin(), score_sums.end()) - score_sums.begin();
        this->star_sequence = sequences[star_pos];
        this->star_sequence_pos = star_pos;
        return this->star_sequence;
    }


public:
    StarAlignment(const std::vector<std::string> & sequences)
    {
        this->sequences = sequences;
    }

    std::string get_star_sequence()
    {
        return this->star_sequence;
    }
    std::vector<std::vector<int>> get_pair_alignment_scores()
    {
        return this->pair_alignment_scores;
    }
    std::vector<std::vector<std::pair<std::string, std::string>>> get_pair_optimal_alignments()
    {
        return this->pair_optimal_alignments;
    }
    int get_star_sequence_pos()
    {
        return this->star_sequence_pos;
    }

    std::vector<std::string> get_optimal_alignment()
    {
        this->find_star_sequence();

        std::vector<std::pair<std::string, std::string>> alignments_with_star = this->pair_optimal_alignments[star_sequence_pos];
        std::vector<std::string> optimal_alignment(sequences.size(), "");
        std::vector<int> pointers(sequences.size(), 0);

        int continous_gaps = 0;
        
        for (int i = 0; i < sequences[star_sequence_pos].size(); i++)
        {
            // La estrella se agrega al inicio
            optimal_alignment[star_sequence_pos] += this->star_sequence[i];
            for (int j = 0; j < sequences.size(); j++)
            {
                // No comparamos con la estrella
                if (j == star_sequence_pos) continue;
                if (pointers[j] >= alignments_with_star[j].second.size()) continue;

                std::string gaps = "";
                bool star_had_gap = false;
                while (alignments_with_star[j].first[pointers[j]] == '_')
                {
                    gaps += '_';
                    continous_gaps++;
                    optimal_alignment[j] += alignments_with_star[j].second[pointers[j]];
                    pointers[j]++;
                    star_had_gap = true;
                }
                if (star_had_gap)
                {
                    for (int p = 0; p < j; p++)
                    {
                        size_t insert_pos = optimal_alignment[p].size() - 1;
                        optimal_alignment[p].insert(insert_pos, gaps);
                    }
                    for (int p = j + 1; p < sequences.size(); p++)
                    {
                        optimal_alignment[p] += gaps;
                    }
                }
                else
                {
                    optimal_alignment[j] += alignments_with_star[j].second[pointers[j]];
                    pointers[j]++;
                }
            }   
        }

        for (int i = 0; i < sequences.size(); i++)
        {
            if (i == star_sequence_pos) continue;
            if (pointers[i] < alignments_with_star[i].second.size())
            {
                optimal_alignment[i] += alignments_with_star[i].second.substr(pointers[i], alignments_with_star[i].second.size() - pointers[i]);
            }
        }
        int max_length = 0;
        for (int i = 0; i < sequences.size(); i++)
        {
            if (optimal_alignment[i].size() > max_length) max_length = optimal_alignment[i].size();
        }

        for (int i = 0; i < sequences.size(); i++)
        {
            while(optimal_alignment[i].size() < max_length) 
                optimal_alignment[i].push_back('_');
        }

        return optimal_alignment;
    }

};

int main()
{   
    std::vector<std::string> sequences = {"ATTGCCATT","ATGGCCATT","ATCCAATTTT","ATCTTCTT","ACTGACC"};
    StarAlignment star(sequences);

    std::ofstream file("result.txt");

    std::vector<std::string> res = star.get_optimal_alignment();
    auto matrix = star.get_pair_alignment_scores();
    auto alignments_with_star =  star.get_pair_optimal_alignments()[star.get_star_sequence_pos()];

    file << "Matriz de Scores:\n\t";
    for (int i = 0; i < sequences.size(); i++)
    {
        file << "S" << i + 1 << "\t";
    }
    file << '\n';

    for (int i = 0; i < matrix.size(); i++)
    {
        file << "S" << i + 1 << "\t";
        int sum = 0;
        for (int j = 0; j < matrix[i].size(); j++)
        {
            if (i == j) file << "-\t";
            else file << matrix[i][j] << '\t';
            sum += matrix[i][j];
        }
        file << "->\t" << sum << '\n';
    }

    file << "Alineamientos de cada secuencia con la estrella:\n";
    for (int i = 0; i < sequences.size(); i++)
    {
        if (i == star.get_star_sequence_pos()) continue;
        file << "S" << star.get_star_sequence_pos() + 1 << ": " << alignments_with_star[i].first << '\n';
        file << "S" << i + 1 << ": " << alignments_with_star[i].second << "\n\n";
    }

    file << "Alineamiento multiple:\n";
    for (int i = 0; i < res.size(); i++)
    {
        file << "S" << i + 1 << ": " << res[i] << '\n';
    }

    file.close();

    return 0;
}
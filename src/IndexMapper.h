#include <vector>
#include <string>
#include <unordered_map>

class IndexMapper {
public:
    static std::unordered_map<std::string, std::vector<int>> mapVector(
        std::vector<std::string>& v
    ) {
        std::unordered_map<std::string, std::vector<int>> map;
        for (int i = 0; i < (int) v.size(); i++) {
            map[v[i]].push_back(i);
        }
        return map;
    }
};

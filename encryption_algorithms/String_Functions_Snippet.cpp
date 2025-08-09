#include <bits/stdc++.h>
using namespace std;



//Declaring Functions
string Normalize(string);
int String_Compare(string, string);
string String_Addition(string, string);
string String_Subtraction(string, string);
string String_Multiplication(string, string);
string String_Divide(string, string);
string String_Modulo(string, string);
string String_GCD(string, string);
string String_Modulo_Exponentiation(string, string, string);
pair<string, pair<string, string>> String_Extended_GCD(string, string);
string String_Modulo_Inverse(string, string);





//Function for making the number usable, removing leading zeros.
string Normalize(string s) {

    //Checking if we have to add a negative sign at the end.
    bool neg = false;
    if (s[0] == '-') {
        neg = true;
        s = s.substr(1);
    }

    //Removing the leading zeros.
    int i = 0;
    while (i < s.size() - 1 && s[i] == '0') i++;
    s = s.substr(i);

    //adding the negative sign back.
    if (neg && s != "0") s = '-' + s;
    return s;
}



//This function returns -1,1 and 0 when s1 is less than, greater than and equal to s2 respectively.
int String_Compare(string s1, string s2) {
    //Remove the leading zeros.
    s1 = Normalize(s1);
    s2 = Normalize(s2);
    
    //Checking if either of them is negative.
    bool neg1 = (s1[0] == '-');
    bool neg2 = (s2[0] == '-');

    //Basic sign check
    if (neg1 && !neg2) return -1; 
    if (!neg1 && neg2) return 1;  
    
    //If both are negative we invert the mod inequality
    if (neg1 && neg2) {
        s1 = s1.substr(1);
        s2 = s2.substr(1);
        if (s1.size() > s2.size()) return -1;
        if (s1.size() < s2.size()) return 1;
        if (s1 > s2) return -1;
        if (s1 < s2) return 1;
        return 0;
    }

    //If both positive we just compare them.
    if (s1.size() > s2.size()) return 1;
    if (s1.size() < s2.size()) return -1;
    if (s1 > s2) return 1;
    if (s1 < s2) return -1;
    return 0;
}



//Function for subtracting the strings.
string String_Subtraction(string s1, string s2) {
    //Making number usable
    s1 = Normalize(s1);
    s2 = Normalize(s2);
    
    //Checking the sign of each of the numbers
    bool neg1 = (s1[0] == '-');
    bool neg2 = (s2[0] == '-');
    
    //Checking whether the difference will be positive or negative and calling functions accordingly.
    if (neg1 && !neg2) return '-' + String_Addition(s1.substr(1), s2);
    if (!neg1 && neg2) return String_Addition(s1, s2.substr(1));
    if (neg1 && neg2)  return String_Subtraction(s2.substr(1), s1.substr(1));
    if (String_Compare(s1, s2) < 0) return '-' + String_Subtraction(s2, s1);

    // Actual subtraction
    reverse(s1.begin(), s1.end());
    reverse(s2.begin(), s2.end());
    
    //Making the length of both numbers same. (Note s1.size() >= s2.size())
    while (s2.size() < s1.size()) s2 += '0';
    
    //Initializing variables
    string result = "";
    int borrow = 0;
    
    //Implimenting the normal subtraction and borrow.
    for (int i = 0; i < s1.size(); i++) {

        //Storing digits.
        int digit1 = s1[i] - '0';
        int digit2 = (i < s2.size()) ? s2[i] - '0' : 0;
        
        //Borrow implimentation.
        int diff = digit1 - digit2 - borrow;
        if (diff < 0) {
            diff += 10;
            borrow = 1;
        } else {
            borrow = 0;
        }
        result += (diff + '0');
    }
    
    // Remove leading zeros
    while (result.size() > 1 && result.back() == '0') {
        result.pop_back();
    }
    
    //Reverting the number back.
    reverse(result.begin(), result.end());
    return Normalize(result);
}



//Function for adding the strings.
string String_Addition(string s1, string s2) {
    //Making the numbers usable.
    s1 = Normalize(s1);
    s2 = Normalize(s2);
    
    //Checking the sign of the numbers.
    bool neg1 = (s1[0] == '-');
    bool neg2 = (s2[0] == '-');
    
    //Checking the sign of the sum and calling functions accordingly.
    if (neg1 && !neg2) return String_Subtraction(s2, s1.substr(1));
    if (!neg1 && neg2) return String_Subtraction(s1, s2.substr(1));
    if (neg1 && neg2) return '-' + String_Addition(s1.substr(1), s2.substr(1));

    // Actual addition algorithm.
    reverse(s1.begin(), s1.end());
    reverse(s2.begin(), s2.end());
    
    //Making the size of both the numbers equal.
    int n = max(s1.size(), s2.size());
    while (s1.size() < n) s1 += '0';
    while (s2.size() < n) s2 += '0';
    
    //Initializing variables
    string result = "";
    int carry = 0;
    
    //Implimentation of digit-wise sum and carry.
    for (int i = 0; i < n; i++) {
        int sum = (s1[i] - '0') + (s2[i] - '0') + carry;
        result += (sum % 10 + '0');
        carry = sum / 10;
    }
    if (carry) result += (carry + '0');
    
    //Reverting the number back.
    reverse(result.begin(), result.end());
    return Normalize(result);
}



//Function of multipying the strings.
string String_Multiplication(string s1, string s2) {
    //Making the numbers useful.
    s1 = Normalize(s1);
    s2 = Normalize(s2);
    
    //If ither one is zero product will be zero.
    if (s1 == "0" || s2 == "0") return "0";
    
    //Checking sign of the numbers
    bool neg1 = (s1[0] == '-');
    bool neg2 = (s2[0] == '-');

    //Checking sign of the product.
    bool result_negative = (neg1 != neg2);
    
    //Taking the absolute value.
    if (neg1) s1 = s1.substr(1);
    if (neg2) s2 = s2.substr(1);
    
    //Initializing Variables.
    int n = s1.size(), m = s2.size();
    vector<int> result(n + m, 0);
    
    //Actual multiplication algorithm
    for (int i = n - 1; i >= 0; i--) {
        for (int j = m - 1; j >= 0; j--) {
            int mul = (s1[i] - '0') * (s2[j] - '0');
            int sum = mul + result[i + j + 1];
            result[i + j + 1] = sum % 10;
            result[i + j] += sum / 10;
        }
    }
    
    //Changing the answer back into a string.
    string ans = "";
    bool leading_zero = true;
    for (int digit : result) {
        if (digit == 0 && leading_zero) continue;
        leading_zero = false;
        ans += (digit + '0');
    }

    if (ans.empty()) ans = "0";
    if (result_negative && ans != "0") ans = '-' + ans;
    return Normalize(ans);
}



//Function for diving the strings. Note that this function does integer division, that is it will not return fractions or decimals.
string String_Divide(string dividend, string divisor) {
    //Making the numbers usable.
    dividend = Normalize(dividend);
    divisor = Normalize(divisor);
    
    //Division by zero is not possible but for simplicity here it just returns 0;
    if (divisor == "0") return "0"; 
    if (dividend == "0") return "0";
    
    //Checking the sign of numbers and the quotient.
    bool neg_dividend = (dividend[0] == '-');
    bool neg_divisor = (divisor[0] == '-');
    bool result_negative = (neg_dividend != neg_divisor);
    
    //Taking the absolute value of the number.
    if (neg_dividend) dividend = dividend.substr(1);
    if (neg_divisor) divisor = divisor.substr(1);
    
    //No fractions in integer division.
    if (String_Compare(dividend, divisor) < 0) return "0";
    
    //Initializing Variables
    string result = "";
    string current = "";
    
    //Actual Division algorithm.
    for (int i = 0; i < dividend.size(); i++) {
        current += dividend[i];
        current = Normalize(current);
        
        int quotient_digit = 0;
        while (String_Compare(current, divisor) >= 0) {
            current = String_Subtraction(current, divisor);
            quotient_digit++;
        }
        
        result += (quotient_digit + '0');
    }
    
    // Remove leading zeros
    int pos = 0;
    while (pos < result.size() - 1 && result[pos] == '0') pos++;
    result = result.substr(pos);
    
    if (result_negative && result != "0") result = '-' + result;
    
    return Normalize(result);
}



//Function for calculating the modulus residue of s1 with respect to s2, which is usually represented as s1%s2.
string String_Modulo(string dividend, string divisor) {
    //Making numbers useable.
    dividend = Normalize(dividend);
    divisor = Normalize(divisor);
    
    //Again division by 0 results in a return value of 0 for simplicity.
    if (divisor == "0") return "0"; 
    if (dividend == "0") return "0";
    
    //Checking sign of numbers
    bool neg_dividend = (dividend[0] == '-');
    bool neg_divisor = (divisor[0] == '-');
    
    //Taking absolute value
    if (neg_dividend) dividend = dividend.substr(1);
    if (neg_divisor) divisor = divisor.substr(1);
    
    //If the mod of s1 is less than s2.
    if (String_Compare(dividend, divisor) < 0) return neg_dividend ? '-' + dividend : dividend;
    
    //Same as the human way of solving the remainder.
    string quotient = String_Divide(dividend, divisor);
    string product = String_Multiplication(quotient, divisor);
    string remainder = String_Subtraction(dividend, product);
    
    //Making remainder positive, just in case.
    if (neg_dividend && remainder != "0") remainder = String_Subtraction(divisor, remainder);
    
    return Normalize(remainder);
}



//Function for finding GCD of the strings
string String_GCD(string s1, string s2) {

    //Normal Euclidian GCD algorithm.
    if (Normalize(s2) == "0") return Normalize(s1);
    return String_GCD(s2, String_Modulo(s1, s2));
}



//Function for calculating modular exponent.{Value = (base)^exponent (Mod mod)}
string String_Modulo_Exponentiation(string base, string exponent, string mod) {

    //Initializing variables.
    string result = "1";
    base = String_Modulo(base, mod);

    //Usual implimentation of the interative Binary Exponentiation.
    while (String_Compare(exponent, "0")) {
        if(exponent == "0") break;
        if ((exponent.back() - '0') % 2 == 1) result = String_Modulo(String_Multiplication(result, base), mod);
        base = String_Modulo(String_Multiplication(base, base), mod);
        exponent = String_Divide(exponent, "2");
    }

    //Returning the value.
    return result;
}



//Function returning the coeffeciets of the Bezout's Lemma. (ax + by = GCD(a,b))
pair<string, pair<string, string>> String_Extended_GCD(string a, string b) {

    //Making numbers usable.
    a = Normalize(a);
    b = Normalize(b);
    
    //Initializing variables.
    string old_r = a, r = b;
    string old_s = "1", s = "0";
    string old_t = "0", t = "1";
    
    //Using the common Bezout's recursive relation to solve the problm iteratively.
    while (r != "0") {
        string quotient = String_Divide(old_r, r);
        
        string new_r = String_Subtraction(old_r, String_Multiplication(quotient, r));
        old_r = r;
        r = new_r;
        
        string new_s = String_Subtraction(old_s, String_Multiplication(quotient, s));
        old_s = s;
        s = new_s;
        
        string new_t = String_Subtraction(old_t, String_Multiplication(quotient, t));
        old_t = t;
        t = new_t;
    }
    
    return make_pair(old_r, make_pair(old_s, old_t));
}



// Calculating the multiplicative inverse of a number in the set of all modular residues of another number.
string String_Modulo_Inverse(string a, string mod) {

    //Making numbers usable.
    a = Normalize(a);
    mod = Normalize(mod);
    
    //Intilializing variables.
    auto result = String_Extended_GCD(a, mod);
    string gcd = result.first;
    string x = result.second.first;
    
    // If they arent coprime, no inverse exists.
    if (gcd != "1") return 0;
    
    // Make sure x is positive and less than mod.
    x = Normalize(x);
    while (x[0] == '-') {
        x = String_Addition(x, mod);
    }
    x = String_Modulo(x, mod);
    
    return Normalize(x);
}



// Function to solve CRT congruencies.
string CRT(vector<string> &Mods, vector<string> &Values) {

    //Total number of congruencies
    long long Length = Mods.size();  
    //Initializing variables
    string M = "1";            
    string CRT_Number = "0";   

    //Finding multiplication of all mods
    for (long long i = 0; i < Length; ++i) {
        M = String_Multiplication(M, Mods[i]);
    }

    //Using the CRT algorithm
    for (long long i = 0; i < Length; ++i) {
        string mi = Mods[i];
        string ai = Values[i];

        string Mi = String_Divide(M, mi); 
        string yi = String_Modulo_Inverse(Mi, mi); 

        
        if (yi == "0") return "-1";  

        
        string term = String_Multiplication(ai, Mi);
        term = String_Multiplication(term, yi);

        
        CRT_Number = String_Addition(CRT_Number, term);
    }

    //Returning the answer mod M
    return String_Modulo(CRT_Number, M);
}



//Function for converting the message into a integer (which is stored as an string.)
string Convert_To_Int(string s) {

    //Initializing variables
    string result = "";

    //Using ASCII values.
    for (char c : s) {
        int ascii = (int)c;
        string padded = string(3 - to_string(ascii).size(), '0') + to_string(ascii);
        result += padded;
    }

    return result;
}



//Function to retrive the message from an integer (which is stored as an string.)
string Convert_To_String(string s) {
    string result = "";
    for (int i = 0; i < s.size(); i += 3) {
        int val = stoi(s.substr(i, 3));
        result += (char)val;
    }
    return result;
}

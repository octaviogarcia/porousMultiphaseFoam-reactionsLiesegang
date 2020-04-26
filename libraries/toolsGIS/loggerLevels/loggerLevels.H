#ifndef LOGGER_LEVELS_H
#define LOGGER_LEVELS_H

extern const int LOG_LEVEL;

template<unsigned int N>
struct LoggerLevel{};

template<unsigned int N,typename T>
const LoggerLevel<N>& operator<<(const LoggerLevel<N>& chain, const T& var){
	if(N <= LOG_LEVEL) Foam::Info << var;
	return chain;
}

#endif
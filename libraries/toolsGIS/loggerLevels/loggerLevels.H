#ifndef LOGGER_LEVELS_H
#define LOGGER_LEVELS_H

extern const int LOG_LEVEL;

template<unsigned int N>
struct loggerLevel{};

template<unsigned int N,typename T>
const loggerLevel<N>& operator<<(const loggerLevel<N>& chain, const T& var){
	if(N <= LOG_LEVEL) Foam::Info << var;
	return chain;
}

#endif
#ifndef ANTOK_EVENT_H
#define ANTOK_EVENT_H

#include<vector>


namespace antok {

	class Function;
	class Initializer;

	class Event {

		friend class antok::Initializer;

	  public:

		static Event* instance();

		bool update();

	  private:

		Event() { };

		static Event* _event;

		std::vector<antok::Function*> _functions;

	};

}

#endif

public:
        void PrintReport( Int thread, float time ) const
        {
            logprint("Thread "+ToString(thread)+": "
                + "\n\t kernel:       \t" + this->  ClassName()
                + "\n\t time elapsed: \t" + ToString(time)
                + "\n"
            );
        }

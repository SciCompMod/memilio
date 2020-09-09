#ifndef EPI_DEREF_ITERATOR_H
#define EPI_DEREF_ITERATOR_H

namespace epi
{

    /**
     * iterator adaptor that makes a range of T* look like a range of T
     * e.g. std::vector<int*> v;
     * then *PointerDereferencingIterator<...>(v.begin()) == **v.begin()
     */
    template<class PtrIter>
    class PointerDereferencingIterator
    {
        private:
        PtrIter m_ptr_iter;

        public:
        using reference = decltype(*std::declval<typename PtrIter::value_type>());
        using value_type = typename std::remove_reference<reference>::type;
        using pointer = value_type*;
        using difference_type = typename PtrIter::difference_type;
        using iterator_category = typename PtrIter::iterator_category;

        PointerDereferencingIterator(PtrIter ptr_iter)
            : m_ptr_iter(ptr_iter)
        {

        }

        reference operator*() const
        {
            return **m_ptr_iter;
        }
        reference operator[](difference_type d) const
        {
            return *m_ptr_iter[d];
        }
        pointer operator->() const
        {
            return m_ptr_iter.operator*().operator*();
        }

        PointerDereferencingIterator& operator++()
        {
            ++m_ptr_iter;
            return *this;
        }
        PointerDereferencingIterator operator++(int)
        {
            return m_ptr_iter++;
        }
        PointerDereferencingIterator& operator+=(difference_type d)
        {
            m_ptr_iter += d;
            return *this;
        }
        PointerDereferencingIterator operator+(difference_type d)
        {
            return m_ptr_iter + d;
        }     
        friend PointerDereferencingIterator operator+(difference_type d, const PointerDereferencingIterator& i)
        {
            return i + d;
        }

        PointerDereferencingIterator& operator--()
        {
            --m_ptr_iter;
            return *this;
        }
        PointerDereferencingIterator operator--(int)
        {
            return m_ptr_iter--;
        }
        PointerDereferencingIterator& operator-=(difference_type d)
        {
            m_ptr_iter -= d;
            return *this;
        }
        PointerDereferencingIterator operator-(difference_type d)
        {
            return m_ptr_iter - d;
        }
        difference_type operator-(PointerDereferencingIterator d)
        {
            return m_ptr_iter - d.m_ptr_iter;
        }

        bool operator<(const PointerDereferencingIterator& other) const
        {
            return m_ptr_iter < other.m_ptr_iter;
        }
        bool operator<=(const PointerDereferencingIterator& other) const
        {
            return m_ptr_iter <= other.m_ptr_iter;
        }
        bool operator>(const PointerDereferencingIterator& other) const
        {
            return m_ptr_iter > other.m_ptr_iter;
        }
        bool operator>=(const PointerDereferencingIterator& other) const
        {
            return m_ptr_iter >= other.m_ptr_iter;
        }
        bool operator!=(const PointerDereferencingIterator& other) const        
        {
            return m_ptr_iter != other.m_ptr_iter;
        }
        bool operator==(const PointerDereferencingIterator& other) const
        {
            return m_ptr_iter == other.m_ptr_iter;
        }
    };
}

#endif
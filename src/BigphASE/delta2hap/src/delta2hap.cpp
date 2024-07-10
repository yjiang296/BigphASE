//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: show-snps.cc
//         Date: 12 / 08 / 2004
//
//        Usage: show-snps [options] <deltafile>
//               Try 'show-snps -h' for more information
//
//  Description: For use in conjunction with the MUMmer package. "show-snps"
//              displays human readable (and machine parse-able) single
//             base-pair polymorphisms, including indels from the .delta output
//            of the "nucmer" program. Outputs SNP positions and relative
//          information to stdout.
//
//------------------------------------------------------------------------------

#include <delta.h>
#include <tigrinc.h>
#include <translate.h>
#include <sw_alignscore.h>
#include <redirect_to_pager.h>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <cstring>
#include <map>
#include <set>
#include <cstdio>
using namespace std;

//=============================================================== Options ====//
string OPT_AlignName;          // delta file name
string OPT_ReferenceName;      // reference sequence file name
string OPT_QueryName;          // query sequence file name
string OPT_prefix = "extract"; // -p option

bool OPT_SortReference = false; // -r option
bool OPT_SortQuery = false;     // -q option
bool OPT_ShowConflict = false;  // -C option
bool OPT_ShowIndels = true;     // -I option
int OPT_window_size = 0;        //-F option
int OPT_max_inclusion = 0;      //-f option

int OPT_padding = 30;
int OPT_max_var_count = 50;
int OPT_Context = 0; // -x option

//============================================================= Constants ====//
const char INDEL_CHAR = '.';
const char SEQEND_CHAR = '-';

struct SNP_R_Sort
{
    bool operator()(const SNP_t *a, const SNP_t *b)
    {
        int i = a->ep->refnode->id->compare(*(b->ep->refnode->id));

        if (i < 0)
            return true;
        else if (i > 0)
            return false;
        else
        {
            if (a->pR < b->pR)
                return true;
            else if (a->pR > b->pR)
                return false;
            else
            {
                int j = a->ep->qrynode->id->compare(*(b->ep->qrynode->id));

                if (j < 0)
                    return true;
                else if (j > 0)
                    return false;
                else
                {
                    if (a->pQ < b->pQ)
                        return true;
                    else
                        return false;
                }
            }
        }
    }
};

struct SNP_Q_Sort
{
    bool operator()(const SNP_t *a, const SNP_t *b)
    {
        int i = a->ep->qrynode->id->compare(*(b->ep->qrynode->id));

        if (i < 0)
            return true;
        else if (i > 0)
            return false;
        else
        {
            if (a->pQ < b->pQ)
                return true;
            else if (a->pQ > b->pQ)
                return false;
            else
            {
                int j = a->ep->refnode->id->compare(*(b->ep->refnode->id));

                if (j < 0)
                    return true;
                else if (j > 0)
                    return false;
                else
                {
                    if (a->pR < b->pR)
                        return true;
                    else
                        return false;
                }
            }
        }
    }
};

// id type chr pos base
struct VAR_t
{
    std::string id;
    std::string type;
    std::string chr;
    long pos;
    std::string base;

    VAR_t()
        : id("un"), type(""), chr(""), pos(0), base("") {}
    void clear()
    {
        id.clear();
        type.clear();
        chr.clear();
        pos = 0;
        base.clear();
    }
};
//========================================================== Fuction Decs ====//
//------------------------------------------------------------------ RevC ----//
inline long RevC(long coord, long len)
{
    return len - coord + 1;
}

//------------------------------------------------------------------ Norm ----//
inline long Norm(long c, long l, int f, AlignmentType_t d)
{
    long retval = (d == PROMER_DATA ? c * 3 - (3 - abs(f)) : c);
    if (f < 0)
        retval = RevC(retval, l);
    return retval;
}

//------------------------------------------------------------------ Swap ----//
inline void Swap(long &a, long &b)
{
    long t = a;
    a = b;
    b = t;
}

//------------------------------------------------------------- CheckSNPs ----//
void CheckSNPs(DeltaGraph_t &graph);

//-------------------------------------------------------------- FindSNPs ----//
void FindSNPs(DeltaGraph_t &graph);

void FilterVars(vector<const VAR_t *> &vars, int max_inclusion, int window_size)
{
    if (max_inclusion < 1)
    {
        throw invalid_argument("max_inclusion must be greater than 0");
    }
    else if (max_inclusion >= window_size)
    {
        return;
    }

    // filter vars
    vector<const VAR_t *> filtered_vars;
    long vars_count = vars.size();

    long cur_var = 0;
    int max_match = 1;
    while (cur_var < vars_count)
    {
        while (cur_var + max_match < vars_count && vars[cur_var + max_match]->pos - vars[cur_var]->pos < window_size && max_match < window_size)
        {
            ++max_match;
        }
        if (max_match > max_inclusion)
        {
            cur_var += max_match;
        }
        else
        {
            VAR_t *varp = new VAR_t;
            *varp = *vars[cur_var];
            filtered_vars.push_back(varp);

            ++cur_var;
        }
        max_match = 1;
    }

    for (vector<const VAR_t *>::iterator it = vars.begin(); it != vars.end(); ++it)
    {
        delete *it;
    }
    vars = filtered_vars;
}

void FindHaplotype(vector<const VAR_t *> &vars, vector<int> &haplotypes, int padding, int max_inclusion)
{
    long vars_count = vars.size();
    long cur_var = 0;
    int max_match = 0;

    haplotypes.clear();

    while (cur_var < vars_count)
    {
        haplotypes.push_back(cur_var);
        while (cur_var + 1 < vars_count && max_match < max_inclusion && vars[cur_var + 1]->pos - vars[cur_var]->pos < padding && vars[cur_var + 1]->pos - vars[cur_var]->pos > 0)
        {
            ++cur_var;
            ++max_match;
        };
        haplotypes.push_back(cur_var);
        max_match = 0;
        ++cur_var;
    }
}
//---------------------------------------------------------- Extract ----//
void Extract(const vector<const SNP_t *> &snps);

//------------------------------------------------------------- ParseArgs ----//
void ParseArgs(int argc, char **argv);

//------------------------------------------------------------- PrintHelp ----//
void PrintHelp(const char *s);

//------------------------------------------------------------ PrintUsage ----//
void PrintUsage(const char *s);

//========================================================= Function Defs ====//
int main(int argc, char **argv)
{
    vector<const SNP_t *> snps;
    DeltaGraph_t graph;

    //-- Command line parsing
    ParseArgs(argc, argv);

    //-- Build the alignment graph from the delta file
    graph.build(OPT_AlignName, true);

    //-- Read sequences
    graph.loadSequences();

    //-- Locate the SNPs
    FindSNPs(graph);

    //-- Check for ambiguous alignment regions
    CheckSNPs(graph);

    //-- Collect and sort the SNPs
    map<string, DeltaNode_t>::iterator mi;
    vector<DeltaEdge_t *>::iterator ei;
    vector<DeltaEdgelet_t *>::iterator li;
    vector<SNP_t *>::iterator si;
    for (mi = graph.refnodes.begin(); mi != graph.refnodes.end(); ++mi)
        for (ei = mi->second.edges.begin(); ei != mi->second.edges.end(); ++ei)
            for (li = (*ei)->edgelets.begin(); li != (*ei)->edgelets.end(); ++li)
                for (si = (*li)->snps.begin(); si != (*li)->snps.end(); ++si)
                    if ((OPT_ShowConflict ||
                         ((*si)->conR == 0 && (*si)->conQ == 0)) &&
                        (OPT_ShowIndels ||
                         ((*si)->cR != INDEL_CHAR && (*si)->cQ != INDEL_CHAR)))
                        snps.push_back(*si);

    if (OPT_SortReference)
        sort(snps.begin(), snps.end(), SNP_R_Sort());
    else
        sort(snps.begin(), snps.end(), SNP_Q_Sort());

    Extract(snps);

    return EXIT_SUCCESS;
}

//------------------------------------------------------------- CheckSNPs ----//
void CheckSNPs(DeltaGraph_t &graph)
{
    map<string, DeltaNode_t>::const_iterator mi;
    vector<DeltaEdge_t *>::const_iterator ei;
    vector<DeltaEdgelet_t *>::iterator eli;
    vector<SNP_t *>::iterator si;
    long i;

    //-- For each reference sequence
    long ref_size = 0;
    long ref_len = 0;
    unsigned char *ref_cov = NULL;
    for (mi = graph.refnodes.begin(); mi != graph.refnodes.end(); ++mi)
    {
        //-- Reset the reference coverage array
        ref_len = (mi->second).len;
        if (ref_len > ref_size)
        {
            ref_cov = (unsigned char *)Safe_realloc(ref_cov, ref_len + 1);
            ref_size = ref_len;
        }
        for (i = 1; i <= ref_len; ++i)
            ref_cov[i] = 0;

        //-- Add to the reference coverage
        for (ei = (mi->second).edges.begin();
             ei != (mi->second).edges.end(); ++ei)
            for (eli = (*ei)->edgelets.begin();
                 eli != (*ei)->edgelets.end(); ++eli)
                for (i = (*eli)->loR; i <= (*eli)->hiR; i++)
                    if (ref_cov[i] < UCHAR_MAX)
                        ref_cov[i]++;

        //-- Set the SNP conflict counter
        for (ei = (mi->second).edges.begin();
             ei != (mi->second).edges.end(); ++ei)
            for (eli = (*ei)->edgelets.begin();
                 eli != (*ei)->edgelets.end(); ++eli)
                for (si = (*eli)->snps.begin(); si != (*eli)->snps.end(); ++si)
                    (*si)->conR = ref_cov[(*si)->pR] - 1;
    }
    free(ref_cov);

    //-- For each query sequence
    long qry_size = 0;
    long qry_len = 0;
    unsigned char *qry_cov = NULL;
    for (mi = graph.qrynodes.begin(); mi != graph.qrynodes.end(); ++mi)
    {
        //-- Reset the query coverage array
        qry_len = (mi->second).len;
        if (qry_len > qry_size)
        {
            qry_cov = (unsigned char *)Safe_realloc(qry_cov, qry_len + 1);
            qry_size = qry_len;
        }
        for (i = 1; i <= qry_len; ++i)
            qry_cov[i] = 0;

        //-- Add to the query coverage
        for (ei = (mi->second).edges.begin();
             ei != (mi->second).edges.end(); ++ei)
            for (eli = (*ei)->edgelets.begin();
                 eli != (*ei)->edgelets.end(); ++eli)
                for (i = (*eli)->loQ; i <= (*eli)->hiQ; i++)
                    if (qry_cov[i] < UCHAR_MAX)
                        qry_cov[i]++;

        //-- Set the SNP conflict counter
        for (ei = (mi->second).edges.begin();
             ei != (mi->second).edges.end(); ++ei)
            for (eli = (*ei)->edgelets.begin();
                 eli != (*ei)->edgelets.end(); ++eli)
                for (si = (*eli)->snps.begin(); si != (*eli)->snps.end(); ++si)
                    (*si)->conQ = qry_cov[(*si)->pQ] - 1;
    }
    free(qry_cov);
}

//-------------------------------------------------------------- FindSNPs ----//
void FindSNPs(DeltaGraph_t &graph)
{
    map<string, DeltaNode_t>::iterator mi;
    vector<DeltaEdge_t *>::iterator ei;
    vector<DeltaEdgelet_t *>::iterator li;
    vector<SNP_t *>::iterator si, psi, nsi;

    //-- For each alignment, identify the SNPs
    for (mi = graph.refnodes.begin(); mi != graph.refnodes.end(); ++mi)
        for (ei = mi->second.edges.begin(); ei != mi->second.edges.end(); ++ei)
        {
            SNP_t *snp;
            int ri, qi;
            char *R[] = {(*ei)->refnode->seq, NULL, NULL, NULL, NULL, NULL, NULL};
            char *Q[] = {(*ei)->qrynode->seq, NULL, NULL, NULL, NULL, NULL, NULL};

            long i;
            long lenR = (*ei)->refnode->len;
            long lenQ = (*ei)->qrynode->len;

            for (li = (*ei)->edgelets.begin(); li != (*ei)->edgelets.end(); ++li)
            {
                long delta;
                int frameR, frameQ, sign; // frame:1 if forward, -1 or 4if reverse ; sign:1 if insertion, -1 if deletion
                long sR, eR, sQ, eQ;
                long rpos, qpos, remain; // remain:length of the remaining sequence
                long rctx, qctx;
                long alenR = lenR;
                long alenQ = lenQ;

                //-- Point the coords the right direction
                frameR = 1;
                if ((*li)->dirR == REVERSE_DIR)
                {
                    sR = RevC((*li)->hiR, lenR);
                    eR = RevC((*li)->loR, lenR);
                    frameR += 3;
                }
                else
                {
                    sR = (*li)->loR;
                    eR = (*li)->hiR;
                }

                frameQ = 1;
                if ((*li)->dirQ == REVERSE_DIR)
                {
                    sQ = RevC((*li)->hiQ, lenQ);
                    eQ = RevC((*li)->loQ, lenQ);
                    frameQ += 3;
                }
                else
                {
                    sQ = (*li)->loQ;
                    eQ = (*li)->hiQ;
                }

                //-- Translate coords to AA if necessary
                if (graph.datatype == PROMER_DATA)
                {
                    alenR /= 3;
                    alenQ /= 3;

                    frameR += (sR + 2) % 3;
                    frameQ += (sQ + 2) % 3;

                    // remeber that eR and eQ point to the last base in the codon
                    sR = (sR + 2) / 3;
                    eR = eR / 3;
                    sQ = (sQ + 2) / 3;
                    eQ = eQ / 3;
                }

                ri = frameR;
                qi = frameQ;

                if (frameR > 3)
                    frameR = -(frameR - 3);
                if (frameQ > 3)
                    frameQ = -(frameQ - 3);

                //-- Generate the sequences if needed
                if (R[ri] == NULL)
                {
                    if (graph.datatype == PROMER_DATA)
                    {
                        R[ri] = (char *)Safe_malloc(alenR + 2);
                        R[ri][0] = '\0';
                        Translate_DNA(R[0], R[ri], ri);
                    }
                    else
                    {
                        R[ri] = (char *)Safe_malloc(alenR + 2);
                        R[ri][0] = '\0';
                        strcpy(R[ri] + 1, R[0] + 1);
                        if ((*li)->dirR == REVERSE_DIR)
                            Reverse_Complement(R[ri], 1, lenR);
                    }
                }
                if (Q[qi] == NULL)
                {
                    if (graph.datatype == PROMER_DATA)
                    {
                        Q[qi] = (char *)Safe_malloc(alenQ + 2);
                        Q[qi][0] = '\0';
                        Translate_DNA(Q[0], Q[qi], qi);
                    }
                    else
                    {
                        Q[qi] = (char *)Safe_malloc(alenQ + 2);
                        Q[qi][0] = '\0';
                        strcpy(Q[qi] + 1, Q[0] + 1);
                        if ((*li)->dirQ == REVERSE_DIR)
                            Reverse_Complement(Q[qi], 1, lenQ);
                    }
                }

                //-- Locate the SNPs
                rpos = sR;
                qpos = sQ;
                remain = eR - sR + 1;

                (*li)->frmR = frameR;
                (*li)->frmQ = frameQ;

                istringstream ss;
                ss.str((*li)->delta);

                while (ss >> delta && delta != 0)
                {
                    sign = delta > 0 ? 1 : -1;
                    delta = labs(delta);

                    //-- For all SNPs before the next indel
                    for (i = 1; i < delta; i++)
                        if (R[ri][rpos++] != Q[qi][qpos++])
                        {
                            if (graph.datatype == NUCMER_DATA &&
                                CompareIUPAC(R[ri][rpos - 1], Q[qi][qpos - 1]))
                                continue;

                            snp = new SNP_t;
                            snp->ep = *ei;
                            snp->lp = *li;
                            snp->pR = Norm(rpos - 1, lenR, frameR, graph.datatype);
                            snp->pQ = Norm(qpos - 1, lenQ, frameQ, graph.datatype);
                            snp->cR = toupper(R[ri][rpos - 1]);
                            snp->cQ = toupper(Q[qi][qpos - 1]);

                            for (rctx = rpos - OPT_Context - 1;
                                 rctx < rpos + OPT_Context; rctx++)
                                if (rctx < 1 || rctx > alenR)
                                    snp->ctxR.push_back(SEQEND_CHAR);
                                else if (rctx == rpos - 1)
                                    snp->ctxR.push_back(snp->cR);
                                else
                                    snp->ctxR.push_back(toupper(R[ri][rctx]));

                            for (qctx = qpos - OPT_Context - 1;
                                 qctx < qpos + OPT_Context; qctx++)
                                if (qctx < 1 || qctx > alenQ)
                                    snp->ctxQ.push_back(SEQEND_CHAR);
                                else if (qctx == qpos - 1)
                                    snp->ctxQ.push_back(snp->cQ);
                                else
                                    snp->ctxQ.push_back(toupper(Q[qi][qctx]));

                            (*li)->snps.push_back(snp);
                        }

                    //-- For the indel
                    snp = new SNP_t;
                    snp->ep = *ei;
                    snp->lp = *li;

                    for (rctx = rpos - OPT_Context; rctx < rpos; rctx++)
                        if (rctx < 1)
                            snp->ctxR.push_back(SEQEND_CHAR);
                        else
                            snp->ctxR.push_back(toupper(R[ri][rctx]));

                    for (qctx = qpos - OPT_Context; qctx < qpos; qctx++)
                        if (qctx < 1)
                            snp->ctxQ.push_back(SEQEND_CHAR);
                        else
                            snp->ctxQ.push_back(toupper(Q[qi][qctx]));

                    if (sign > 0)
                    {
                        snp->pR = Norm(rpos, lenR, frameR, graph.datatype);
                        if (frameQ > 0)
                            snp->pQ = Norm(qpos - 1, lenQ, frameQ, graph.datatype);
                        else
                            snp->pQ = Norm(qpos, lenQ, frameQ, graph.datatype);

                        snp->cR = toupper(R[ri][rpos++]);
                        snp->cQ = INDEL_CHAR;

                        remain -= i;
                        rctx++;
                    }
                    else
                    {
                        snp->pQ = Norm(qpos, lenQ, frameQ, graph.datatype);
                        if (frameR > 0)
                            snp->pR = Norm(rpos - 1, lenR, frameR, graph.datatype);
                        else
                            snp->pR = Norm(rpos, lenR, frameR, graph.datatype);

                        snp->cR = INDEL_CHAR;
                        snp->cQ = toupper(Q[qi][qpos++]);

                        remain -= i - 1;
                        qctx++;
                    }

                    snp->ctxR.push_back(snp->cR);
                    for (; rctx < rpos + OPT_Context; rctx++)
                        if (rctx > alenR)
                            snp->ctxR.push_back(SEQEND_CHAR);
                        else
                            snp->ctxR.push_back(toupper(R[ri][rctx]));

                    snp->ctxQ.push_back(snp->cQ);
                    for (; qctx < qpos + OPT_Context; qctx++)
                        if (qctx > alenQ)
                            snp->ctxQ.push_back(SEQEND_CHAR);
                        else
                            snp->ctxQ.push_back(toupper(Q[qi][qctx]));

                    (*li)->snps.push_back(snp);
                }

                //-- For all SNPs after the final indel
                for (i = 0; i < remain; i++)
                    if (R[ri][rpos++] != Q[qi][qpos++])
                    {
                        if (graph.datatype == NUCMER_DATA &&
                            CompareIUPAC(R[ri][rpos - 1], Q[qi][qpos - 1]))
                            continue;

                        snp = new SNP_t;
                        snp->ep = *ei;
                        snp->lp = *li;
                        snp->pR = Norm(rpos - 1, lenR, frameR, graph.datatype);
                        snp->pQ = Norm(qpos - 1, lenQ, frameQ, graph.datatype);
                        snp->cR = toupper(R[ri][rpos - 1]);
                        snp->cQ = toupper(Q[qi][qpos - 1]);

                        for (rctx = rpos - OPT_Context - 1;
                             rctx < rpos + OPT_Context; rctx++)
                            if (rctx < 1 || rctx > alenR)
                                snp->ctxR.push_back(SEQEND_CHAR);
                            else if (rctx == rpos - 1)
                                snp->ctxR.push_back(snp->cR);
                            else
                                snp->ctxR.push_back(toupper(R[ri][rctx]));

                        for (qctx = qpos - OPT_Context - 1;
                             qctx < qpos + OPT_Context; qctx++)
                            if (qctx < 1 || qctx > alenQ)
                                snp->ctxQ.push_back(SEQEND_CHAR);
                            else if (qctx == qpos - 1)
                                snp->ctxQ.push_back(snp->cQ);
                            else
                                snp->ctxQ.push_back(toupper(Q[qi][qctx]));

                        (*li)->snps.push_back(snp);
                    }

                //-- Sort SNPs and calculate distances
                if (OPT_SortReference)
                {
                    sort((*li)->snps.begin(), (*li)->snps.end(), SNP_R_Sort());

                    for (si = (*li)->snps.begin(); si != (*li)->snps.end(); ++si)
                    {
                        psi = si - 1;
                        nsi = si + 1;

                        (*si)->buff = 1 +
                                      ((*si)->pR - (*li)->loR < (*li)->hiR - (*si)->pR ? (*si)->pR - (*li)->loR : (*li)->hiR - (*si)->pR);

                        if (psi >= (*li)->snps.begin() &&
                            (*si)->pR - (*psi)->pR < (*si)->buff)
                            (*si)->buff = (*si)->pR - (*psi)->pR;

                        if (nsi < (*li)->snps.end() &&
                            (*nsi)->pR - (*si)->pR < (*si)->buff)
                            (*si)->buff = (*nsi)->pR - (*si)->pR;
                    }
                }
                else
                {
                    sort((*li)->snps.begin(), (*li)->snps.end(), SNP_Q_Sort());

                    for (si = (*li)->snps.begin(); si != (*li)->snps.end(); ++si)
                    {
                        psi = si - 1;
                        nsi = si + 1;

                        (*si)->buff = 1 +
                                      ((*si)->pQ - (*li)->loQ < (*li)->hiQ - (*si)->pQ ? (*si)->pQ - (*li)->loQ : (*li)->hiQ - (*si)->pQ);

                        if (psi >= (*li)->snps.begin() &&
                            (*si)->pQ - (*psi)->pQ < (*si)->buff)
                            (*si)->buff = (*si)->pQ - (*psi)->pQ;

                        if (nsi < (*li)->snps.end() &&
                            (*nsi)->pQ - (*si)->pQ < (*si)->buff)
                            (*si)->buff = (*nsi)->pQ - (*si)->pQ;
                    }
                }
            }

            //-- Clear up the seq
            for (i = 1; i <= 6; i++)
            {
                free(R[i]);
                free(Q[i]);
            }
        }
}

void Extract(const std::vector<const SNP_t *> &snps)
{
    vector<const SNP_t *>::const_iterator si = snps.begin();
    vector<const VAR_t *> vars; // var(variation) is a snp or an indel
    vector<int> haplotypes;
    // catch all vars
    while (si != snps.end())
    {
        VAR_t *var = new VAR_t;
        var->chr = *((*si)->ep->refnode->id);

        if ((*si)->cR != '.' && (*si)->cQ != '.' && (*si)->cR != (*si)->cQ)
        {
            var->type = "single";
            var->pos = (*si)->pR;
            var->base = (*si)->cQ;
            ++si;
        }
        else if ((*si)->cR == '.' && (*si)->cQ != '.')
        {
            var->type = "insertion";
            var->pos = (*si)->pR;
            int i = 0;
            do
            {
                var->base += (*(si + i))->cQ;
                ++i;
            } while ((si + i) != snps.end() && (*(si + i))->pQ == (*si)->pQ + i && (*(si + i))->cR == '.');
            si += i;
        }
        else if ((*si)->cR != '.' && (*si)->cQ == '.')
        {
            var->type = "deletion";
            var->pos = (*si)->pR;
            int i = 0;
            do
            {
                ++i;
            } while ((si + i) != snps.end() && (*(si + i))->pR == (*si)->pR + i && (*(si + i))->cQ == '.');
            var->base = std::to_string(i);
            si += i;
        }
        else
        {
            std::cerr << "error" << std::endl;
            assert(false);
        }
        vars.push_back(var);
        // std::cout << var.id << "\t" << var.type << "\t" << var.chr << "\t" << var.pos << "\t" << var.base << std::endl;
    }
    if (OPT_max_inclusion && OPT_window_size)
    {
        try
        {
            FilterVars(vars, OPT_max_inclusion, OPT_window_size);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    FindHaplotype(vars, haplotypes, OPT_padding, OPT_max_var_count); // padding 30, max_inclusion 50

    // print vars to a file
    ofstream sofs(OPT_prefix + ".snp");
    for (long i = 0; i < vars.size(); ++i)
    {
        sofs << vars[i]->id + to_string(i) << "\t" << vars[i]->type << "\t" << vars[i]->chr << "\t" << vars[i]->pos << "\t" << vars[i]->base << std::endl;
    }
    // print haplotypes to a file
    ofstream hofs(OPT_prefix + ".haplotype");
    long ht = 0;
    string vars_id;
    for (long i = 0; i < haplotypes.size(); i += 2)
    {
        long var_start = haplotypes[i];
        long var_end = haplotypes[i + 1];
        vars_id.clear();
        vars_id = "un" + to_string(var_start);
        long var_curr = var_start;
        while (++var_curr <= var_end)
        {
            vars_id += ",un" + to_string(var_curr);
        }
        // id chr start end vars
        hofs << "ht" + to_string(ht) << "\t" << vars[var_start]->chr << "\t" << vars[var_start]->pos << "\t" << vars[var_end]->pos << "\t" << vars_id << std::endl;
        ++ht;
    }

    sofs.close();
    hofs.close();
}

//------------------------------------------------------------- ParseArgs ----//
void ParseArgs(int argc, char **argv)
{
    int ch, errflg = 0;
    optarg = NULL;

    while (!errflg &&
           ((ch = getopt(argc, argv, "ChIl:L:p:qrf:F:x:")) != EOF))
        switch (ch)
        {
        case 'C':
            OPT_ShowConflict = true;
            break;

        case 'h':
            PrintHelp(argv[0]);
            exit(EXIT_SUCCESS);
            break;

        case 'I':
            OPT_ShowIndels = false;
            break;

        case 'l':
            OPT_padding = atoi(optarg);
            break;

        case 'L':
            OPT_max_var_count = atoi(optarg);
            break;

        case 'p':
            OPT_prefix = optarg;
            break;

        case 'q':
            OPT_SortQuery = true;
            break;

        case 'r':
            OPT_SortReference = true;
            break;

        case 'f':
            OPT_max_inclusion = atoi(optarg);
            break;
        case 'F':
            OPT_window_size = atoi(optarg);
            break;
        case 'x':
            OPT_Context = atoi(optarg);
            break;

        default:
            errflg++;
        }

    if (OPT_Context < 0)
    {
        cerr << "ERROR: SNP context must be a positive int\n";
        errflg++;
    }

    if (OPT_SortReference && OPT_SortQuery)
        cerr << "WARNING: both -r and -q were passed, -q ignored\n";

    if (!OPT_SortReference && !OPT_SortQuery)
        OPT_SortReference = true;

    if (errflg > 0 || optind != argc - 1)
    {
        PrintUsage(argv[0]);
        cerr << "Try '" << argv[0] << " -h' for more information.\n";
        exit(EXIT_FAILURE);
    }

    OPT_AlignName = argv[optind++];
}

//------------------------------------------------------------- PrintHelp ----//
void PrintHelp(const char *s)
{
    PrintUsage(s);
    cerr
        << "-C            report SNPs from alignments with an ambiguous\n"
        << "              mapping\n"
        << "-h            Display help information\n"
        << "-I            Do not report indels\n"
        << "-l int           maximum distance(bp),default(30) between variations(snp/indel) in a sigle haplotype\n"
        << "-L int           maxmium count,default(50) of variation(snp/indel) a sigle haplotype can have \n"
        << "-p            the output file name(prefix.snp & prefix.haplotype)\n"
        << "-q            Sort output lines by query IDs and SNP positions\n"
        << "-r            Sort output lines by reference IDs and SNP positions\n"
        << "              'show-coords' lines to stdin\n"
        << "-f int           filter varations gereter tahn f (maximum_inclusion) in F (window_size bp) \n"
        << "-F int           filter varations gereter tahn f (maximum_inclusion) in F (window_size bp) \n"
        << "-x int        Include x characters of surrounding SNP context in the\n"
        << "              output, default "
        << OPT_Context << endl
        << endl;

    cerr
        << "  Input is the .delta output of the nucmer program\n"
        << "passed on the command line.\n"
        << "  Output is to two files named prefix.snp & prefix.haplotype \n"
        << endl;

    return;
}

//------------------------------------------------------------ PrintUsage ----//
void PrintUsage(const char *s)
{
    cerr
        << "\nUSAGE: " << s << "  [options]  <deltafile>\n\n";
    return;
}

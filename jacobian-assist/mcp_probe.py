"""Actual stdio MCP initialize/find/run smoke test using pinned Jacobian."""
from pathlib import Path
import argparse,asyncio,json,sys

from mcp import ClientSession,StdioServerParameters
from mcp.client.stdio import stdio_client

HERE=Path(__file__).resolve().parent


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def dump(result):return result.model_dump(mode='json',by_alias=True)


async def run():
    native=json.loads((HERE/'native_report.json').read_bytes())
    case=native['native_Gram_inertia_crosschecks'][0]
    params=StdioServerParameters(command=sys.executable,args=['-m','jacobian.mcp.stdio_cli'],
                                env={'PYTHONUTF8':'1'})
    async with stdio_client(params) as (read,write):
        async with ClientSession(read,write,read_timeout_seconds=60) as session:
            init=await session.initialize()
            need(init.server_info.name=='jacobian' and init.server_info.version=='0.23.0','pinned MCP server')
            tools=await session.list_tools()
            names=[t.name for t in tools.tools]
            need(set(names)=={'math.find','math.run'},'actual mathematical discovery/execution tools')
            contract=await session.call_tool('math.find',{'operation_id':'matrix.inertia.compute'})
            need(not contract.is_error,'contract inspection over MCP')
            result=await session.call_tool('math.run',{'operation_id':'matrix.inertia.compute','payload':case['input']})
            need(not result.is_error,'actual stored Gram block executed over MCP')
            actual=json.loads(next(b.text for b in result.content if b.type=='text'))['output']
            expected=case['output']
            need(all(actual[k]==expected[k] for k in ('n_positive','n_negative','n_zero','definiteness')),
                 'MCP exact inertia agrees with the independently checked native result')
    return {'passed':True,'jacobian_version':'0.23.0','transport':'stdio',
            'server':dump(init)['serverInfo'],'tool_names':names,
            'inspected_operation':'matrix.inertia.compute',
            'operation_input':case['input'],'operation_output':dump(result),
            'MCP_inertia_matches_independently_checked_native_result':True,
            'scope':'real local MCP handshake, contract inspection and operation execution; no global bound claim'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--record',action='store_true')
    args=parser.parse_args();report=asyncio.run(run())
    if args.record:(HERE/'mcp_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
